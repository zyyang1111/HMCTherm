#include "Thermal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> //stringstream
#include <stdlib.h> // getenv()
#include <math.h>




using namespace std; 
using namespace CasHMC; 

extern "C" double ***steady_thermal_solver(double ***powerM, double W, double Lc, int numP, int dimX, int dimZ, double **Midx, int count);
extern "C" double *transient_thermal_solver(double ***powerM, double W, double L, int numP, int dimX, int dimZ, double **Midx, int MidxSize, double *Cap, int CapSize, double time, int iter, double *T_trans);
extern "C" double **calculate_Midx_array(double W, double Lc, int numP, int dimX, int dimZ, int* MidxSize);
extern "C" double *calculate_Cap_array(double W, double Lc, int numP, int dimX, int dimZ, int* CapSize); 
extern "C" double *initialize_Temperature(double W, double Lc, int numP, int dimX, int dimZ); 
extern "C" double get_maxT(double *Tc, int Tsize);

extern string logicPFileName; // the name of the file given the power profiles for the logic layer(s)
extern string resultdir; // the directory name for storing the result
extern long PowerEpoch; 
// used for control the start point when the simulation is restarted
// clk_cycle_dist is the started clock cycle 
// cont_bool = 0 indicates starting from the very beginning (i.e. clk_cycle_dist = 0)
// cont_bool = 1 indicates starting from clk_cycle_dist --> all the previous data will be loaded
extern uint64_t clk_cycle_dist; 
extern int cont_bool; 
extern int num_refresh_save; 

ThermalCalculator::ThermalCalculator(bool withLogic_):
	totalEnergy(0.0),
	sampleEnergy(0.0),
	pe_crit(false),
	sample_id(0),
	withLogic(withLogic_),
	RefreshCont(RFControl())
	{
		num_refresh = num_refresh_save; 
		// if start from the scratch, num_refresh_save = 0
		power_epoch = PowerEpoch; 
		std::cout << "enter the assignment method\n";
		std::cout << "ARCH_SCHEME = " << ARCH_SCHEME << std::endl;

		NUM_GRIDS_X = NUM_ROWS / MAT_X; 
		NUM_GRIDS_Y = NUM_COLS / MAT_Y; 

		if (ARCH_SCHEME == 0)
		{
			vault_x = square_array(NUM_VAULTS); 
			bank_y = square_array(NUM_BANKS); 

			vault_y = NUM_VAULTS / vault_x; 
			bank_x = NUM_BANKS / bank_y;
		}
		else if (ARCH_SCHEME == 1)
		{
			int num_bank_per_layer = NUM_BANKS / NUM_LAYERS; 
			bank_y = num_bank_per_layer; 
			bank_x = 1; 
			vault_x = 1; 
			vault_y = NUM_VAULTS; 
			x = vault_x * bank_x * NUM_GRIDS_X; 
			y = vault_y * bank_y * NUM_GRIDS_Y;
			double asr = max(x, y) / min(x,y); 
			int vault_x_r = vault_x; 
			while (true){
				vault_x ++; 
				vault_y = NUM_VAULTS/vault_x; 
				if (vault_x * vault_y != NUM_VAULTS)
					continue; 
				x = vault_x * bank_x * NUM_GRIDS_X; 
				y = vault_y * bank_y * NUM_GRIDS_Y;
				double asr_n = max(x, y) / min(x,y);
				if (asr_n >= asr)
					break; 
				vault_x_r = vault_x; 
				asr = asr_n;  
			}
			vault_x = vault_x_r; 
			vault_y = NUM_VAULTS / vault_x; 

			//std::cout << "vault_x = " << vault_x << "; vault_y = " << vault_y << "\n"; 
			//std::cout << "bank_x = " << bank_x << "; bank_y = " << bank_y << "\n";
			//std::cout << "NUM_GRIDS_X = " << NUM_GRIDS_X << "; NUM_GRIDS_Y = " << NUM_GRIDS_Y << std::endl;
		}
 		else if (ARCH_SCHEME == 2)
		{
			bank_x = 1; 
			bank_y = 1; 
			vault_x = square_array(NUM_VAULTS); 
			vault_y = NUM_VAULTS / vault_x;
		}

		x = vault_x * bank_x * NUM_GRIDS_X; 
		y = vault_y * bank_y * NUM_GRIDS_Y; 
		z = NUM_LAYERS;

		// initialize the accumulative power maps 
		accu_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
		accu_Pmap_wLogic = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z+2, 0))); // memory logic layer + processor layer
		cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
		cur_Pmap_wLogic = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z+2, 0)));

		vault_usage_single = vector<int> (NUM_VAULTS, 0);
		vault_usage_multi = vector<int> (NUM_VAULTS, 0);

		std::cout << "(x, y, z) = " << "( " << x << ", " << y << ", " << z << " )" << std::endl; 

/*		for (int l = 0; l < z; l ++)
		{
			std::cout << "layer " << l << std::endl; 
			for (int i = 0; i < x; i ++)
			{
				for (int j = 0; j < y; j ++)
				{
					std::cout << accu_Pmap[i][j][l] << " "; 
				}
				std::cout << std::endl;
			}
		}
*/
		std::cout << "vault_x = " << vault_x << "; vault_y = " << vault_y << std::endl;
		std::cout << "bank_x = " << bank_x << "; bank_y = " << bank_y << std::endl;

		///////// calculate Midx and MidxSize for temperature ///////////////
		calcMidx();
		///////// read power of the logic layer from the file //////////
		ReadlogicP();
		///////// calculate Midx and MidxSize for PDN ///////////////
		// calcPDNMidx();
		//////// Initialize the transient PDN variables /////////////
		// IniTransPDN();

		/* define the file name */
		power_trace_str = resultdir + "power_trace.csv"; 
		temp_trace_str = resultdir + "temperature_trace.csv"; 
		avg_power_str = resultdir + "Average_Power_Profile.csv";
		final_temp_str = resultdir + "static_temperature.csv"; 
		debug_power_resize_str = resultdir + "Debug_power_profile_resize.csv";
		debug_power_str = resultdir + "Debug_power_profile.csv"; 
		dump_RT_str = resultdir + "RetTCountDown_file.txt";
		dump_curCyc_str = resultdir + "currentClockCycle_file.txt";
		dump_Ttrans_str = resultdir + "Ttrans_file.txt"; 
		dump_accuP_str = resultdir + "accuP_file.txt";
		dump_curP_str = resultdir + "curP_file.txt";


		/* print the header for csv files */
		//std::ofstream power_file; 
		//std::ofstream temp_file; 
		//power_file.open(power_trace_str.c_str()); power_file << "S_id,layer,x,y,power\n"; power_file.close();
		//temp_file.open(temp_trace_str.c_str()); temp_file << "S_id,layer,x,y,temperature\n"; temp_file.close();

		// if cont_bool == 1, reload the PTdata: T_trans, accu_Pmap, cur_Pmap
		if (cont_bool)
			Reload_PTdata(); 
		else{
			/* print the header for csv files */
			std::ofstream power_file; 
			std::ofstream temp_file; 
			power_file.open(power_trace_str.c_str()); power_file << "S_id,layer,x,y,power\n"; power_file.close();
			temp_file.open(temp_trace_str.c_str()); temp_file << "S_id,layer,x,y,temperature\n"; temp_file.close();
		}

		t = clock();
	}

ThermalCalculator::~ThermalCalculator()
{
	std::cout << "delete ThermalCalculator\n";
	Dump_PTdata(); 
	/* free the space of T */
    for (size_t i = 0; i < x; i++)
    {
        for (size_t j = 0; j < y; j++)
        {
            free(T_final[i][j]);
        }
        free(T_final[i]); 
    }
    free(T_final);

    for (int k = 0; k < MidxSize; k ++)
    {
        free(Midx[k]);
    }
    free(Midx); 
    

    /*for (size_t i = 0; i < PDN_x; i++)
    {
        for (size_t j = 0; j < PDN_y; j++)
        {
            free(V_final[i][j]);
        }
        free(V_final[i]); 
    }
    free(V_final);*/

    for (int k = 0; k < PDNMidxSize; k ++)
    {
        free(PDNMidx[k]);
    }

    free(T_trans); 
    free(Cap);
}

void ThermalCalculator::Dump_PTdata()
{
	std::ofstream Trans_file, accuP_file, curP_file;  
	Trans_file.open(dump_Ttrans_str.c_str());
	accuP_file.open(dump_accuP_str.c_str());
	curP_file.open(dump_curP_str.c_str()); 

	int numP = ( withLogic ? z+2 : z); 
	
	for (int i = 0; i < x * y * (numP*3+1); i++)
		Trans_file << T_trans[i] << " "; 
	Trans_file.close(); 

	for (int i = 0; i < x; i ++){
		for (int j = 0; j < y; j ++){
			for (int k = 0; k < z; k ++){
				accuP_file << accu_Pmap[i][j][k] << " "; 
				curP_file << cur_Pmap[i][j][k] << " "; 
			}
		}
	}
	accuP_file.close(); 
	curP_file.close();
}

void ThermalCalculator::Reload_PTdata()
{
	// The file is the same as the dump file in Dump_PTdata
	std::ifstream Trans_file, accuP_file, curP_file; 
	Trans_file.open(dump_Ttrans_str.c_str());
	accuP_file.open(dump_accuP_str.c_str());
	curP_file.open(dump_curP_str.c_str()); 

	int numP = ( withLogic ? z+2 : z);

	for (int i = 0; i < x * y * (numP*3+1); i ++)
		Trans_file >> T_trans[i]; 
	Trans_file.close(); 

	for (int i = 0; i < x; i ++){
		for (int j = 0; j < y; j ++){
			for (int k = 0; k < z; k ++){
				accuP_file >> accu_Pmap[i][j][k]; 
				curP_file >> cur_Pmap[i][j][k]; 
			}
		}
	}
	accuP_file.close(); 
	curP_file.close();
}


int ThermalCalculator::square_array(int total_grids_)
{
	int x, y, x_re = 1; 
	for (x = 1; x <= sqrt(total_grids_); x ++)
	{
		y = total_grids_ / x; 
		if (x * y == total_grids_)
			x_re = x; 
	}
	return x_re; 
}

void ThermalCalculator::addPower_refresh(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, uint64_t cur_cycle)
{
	if (cur_cycle <= clk_cycle_dist){
		if (cur_cycle > (sample_id+1) * power_epoch)
			sample_id = sample_id + 1;
		return; 
	}
	cout << "addPower_refresh\n";
	num_refresh ++;
	if (cur_cycle > (sample_id+1) * power_epoch)
	{
		save_sampleP(cur_cycle, sample_id); 
		cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
		sampleEnergy = 0; 
		sample_id = sample_id + 1; 
	}
	totalEnergy += energy_t_ * Vdd / 1000.0;
    sampleEnergy += energy_t_ * Vdd / 1000.0;
    
    vault_usage_multi[vault_id_] ++; 
    double energy; 
    int ref_layer = 0, ref_row_phy = 0, ref_col_phy = 0; 
    if (ARCH_SCHEME == 1)
    {
    	energy = energy_t_ / (REFRESH_ROWNUM); 
    	for (int i = 0; i < REFRESH_ROWNUM; i ++)
    	{
    		mapPhysicalLocation(vault_id_, bank_id_, row_id_+i, col_id_, &ref_layer, &ref_row_phy, &ref_col_phy);
    		accu_Pmap[ref_row_phy][ref_col_phy][ref_layer] += energy * Vdd / 1000.0;
    		cur_Pmap[ref_row_phy][ref_col_phy][ref_layer] += energy * Vdd / 1000.0;
    	}
    }

}


void ThermalCalculator::addPower(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, bool single_bank, uint64_t cur_cycle)
{
	//cout << "addPower\n";
	//cout << "cur_cycle = " << cur_cycle << endl;

	//std::cout << "energy = " << energy_t_ << std::endl; 
	//std::cout << "(vault, bank, row, col) = " << "( " << vault_id_ << ", " << bank_id_ << ", " << row_id_ << ", " << col_id_ << " )" << std::endl;
	//std::cout << "single_bank is " << single_bank << std::endl;
    if (cur_cycle <= clk_cycle_dist){
		if (cur_cycle > (sample_id+1) * power_epoch)
			sample_id = sample_id + 1;
		return; 
	}

	////// determine whether the sampling period ends //////////////
	if (cur_cycle > (sample_id+1) * power_epoch)
	{
		save_sampleP(cur_cycle, sample_id); 
		cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
		sampleEnergy = 0; 
		sample_id = sample_id + 1; 
	}
/*
	if (cur_cycle % power_epoch == 0)
		pe_crit = true; 
	else
	{
		if (pe_crit)
		{
			// save the sampling power 
			save_sampleP(cur_cycle); 
			cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
			sampleEnergy = 0;
			pe_crit = false; 
		}
	}*/
	////////////////////////////////////////////////////////////////


    totalEnergy += energy_t_ * Vdd / 1000.0;
    sampleEnergy += energy_t_ * Vdd / 1000.0;
    
	if (single_bank)
	{
		vault_usage_single[vault_id_] ++;
		int layer = 0, row_phy = 0, col_phy = 0; 
		mapPhysicalLocation(vault_id_, bank_id_, row_id_, col_id_, &layer, &row_phy, &col_phy);
		accu_Pmap[row_phy][col_phy][layer] += energy_t_ * Vdd / 1000.0; 
		cur_Pmap[row_phy][col_phy][layer] += energy_t_ * Vdd / 1000.0; 

	}
	else
	{
		vault_usage_multi[vault_id_] ++;
		double energy; 
		int base_layer = 0, base_row_phy = 0, base_col_phy = 0; 
		mapPhysicalLocation(vault_id_, bank_id_, row_id_, col_id_, &base_layer, &base_row_phy, &base_col_phy);
		//std::cout << "(layer, row, col) = " << "( " << base_layer << ", " << base_row_phy << ", " << base_col_phy << " )" << std::endl; 
		
        if (ARCH_SCHEME == 1)
        {
        	energy = energy_t_ / (NUM_GRIDS_X * NUM_GRIDS_Y); 
			for (int i = base_row_phy; i < base_row_phy + NUM_GRIDS_X; i ++)
			{
				for (int j = base_col_phy; j < base_col_phy + NUM_GRIDS_Y; j ++)
				{
					accu_Pmap[i][j][base_layer] += energy * Vdd / 1000.0; 
					cur_Pmap[i][j][base_layer] += energy * Vdd / 1000.0; 
					//std::cout << "(" << i << "," << j << "," << base_layer << ") ADDS " << energy * Vdd / 1000.0 << " NOW IS " << accu_Pmap[i][j][base_layer] << std::endl << std::endl;
				}
			}
        }
        else if (ARCH_SCHEME == 0)
        {
        	energy = energy_t_ / (NUM_GRIDS_X * NUM_GRIDS_Y * NUM_LAYERS); 
			for (int l = 0; l < NUM_LAYERS; l ++){
				for (int i = base_row_phy; i < base_row_phy + NUM_GRIDS_X; i ++){
					for (int j = base_col_phy; j < base_col_phy + NUM_GRIDS_Y; j ++){
						accu_Pmap[i][j][l] += energy * Vdd / 1000.0; 
						cur_Pmap[i][j][l] += energy * Vdd / 1000.0; 
					}
				}
			}
        }	
        else if (ARCH_SCHEME == 2)
        {
        	energy = energy_t_ / NUM_LAYERS; 
        	for (int l = 0; l < NUM_LAYERS; l ++){
        		accu_Pmap[base_row_phy][base_col_phy][l] += energy * Vdd / 1000.0;
        		cur_Pmap[base_row_phy][base_col_phy][l] += energy * Vdd / 1000.0;
        	}
        }	

	}

}

void ThermalCalculator::addIOPower(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, uint64_t cur_cycle)
{
	if (cur_cycle <= clk_cycle_dist){
		if (cur_cycle > (sample_id+1) * power_epoch)
			sample_id = sample_id + 1;
		return; 
	}
	////// determine whether the sampling period ends //////////////
	if (cur_cycle > (sample_id+1) * power_epoch)
	{
		save_sampleP(cur_cycle, sample_id); 
		cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
		sampleEnergy = 0; 
		sample_id = sample_id + 1; 
	}
	/*
	if (cur_cycle % power_epoch == 0)
		pe_crit = true; 
	else
	{
		if (pe_crit)
		{
			// save the sampling power 
			save_sampleP(cur_cycle); 
			cur_Pmap = vector<vector<vector<double> > > (x, vector<vector<double> > (y, vector<double> (z, 0)));
			sampleEnergy = 0;
			pe_crit = false; 
		}
	}*/

	int NGperLperV = bank_x * bank_y * NUM_GRIDS_X * NUM_GRIDS_Y;
	double energy = energy_t_ / NGperLperV; // spread the power to all the grids on one layer
	int layer = 0, row_phy = 0, col_phy = 0; 
	mapPhysicalLocation(vault_id_, bank_id_, row_id_, col_id_, &layer, &row_phy, &col_phy);
	int vault_id_x = vault_id_ / vault_y; 
	int vault_id_y = vault_id_ % vault_y; 

	totalEnergy += energy_t_ * Vdd / 1000.0; 
	IOEnergy += energy_t_ * Vdd / 1000.0;

	// just need the layer of the visit 
	for (int i = vault_id_x * bank_x * NUM_GRIDS_X; i < (vault_id_x + 1) * bank_x * NUM_GRIDS_X; i ++){
		for (int j = vault_id_y * bank_y * NUM_GRIDS_Y; j < (vault_id_y + 1) * bank_y * NUM_GRIDS_Y; j ++){
			accu_Pmap[i][j][layer] += energy * Vdd / 1000.0;
			cur_Pmap[i][j][layer] += energy * Vdd / 1000.0;
		}
	}
}


void ThermalCalculator::rev_mapPhysicalLocation(int *vault_id_, int *bank_id_, int *row_s, int *row_e, int layer, int row, int col)
{
	// currently only support ARCH_SCHEME 1 
	if (ARCH_SCHEME == 1)
	{
		int grid_step = NUM_ROWS / (NUM_GRIDS_X * NUM_GRIDS_Y);
		int num_bank_per_layer = NUM_BANKS / NUM_LAYERS;  
		int bx = row / NUM_GRIDS_X; 
		int by = col / NUM_GRIDS_Y; 
		int vx = bx / bank_x; 
		int vy = by / bank_y; 
		// get the vault id
		*vault_id_ = vx * vault_y + vy; 

		// get the bank id
		int bank_same_layer = (bx % bank_x)* bank_y + (by % bank_y); 
		*bank_id_ = layer * num_bank_per_layer + bank_same_layer; 

		// get the row id
		int grid_id = (row % NUM_GRIDS_X) * NUM_GRIDS_Y + (col % NUM_GRIDS_Y); 
		*row_s = grid_id * grid_step; 
		*row_e = (grid_id+1) * grid_step; 
	}
	else
		cout << "Currently Only Support ARCH_SCHEME = 1\n";
}



void ThermalCalculator::mapPhysicalLocation(unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, int *layer, int *row, int *col)
{

	int vault_id_x = vault_id_ / vault_y; 
	int vault_id_y = vault_id_ % vault_y; 

	if (ARCH_SCHEME == 0)
	{
		// layer # is determined by the index of rows and columns 
		// each bank is divided into NUM_GRIDS_X * NUM_GRIDS_Y * NUM_LAYERS thermal grids 
		*layer = row_id_ / (NUM_ROWS / NUM_LAYERS); 
		int grid_id = (row_id_ % (NUM_ROWS / NUM_LAYERS)) / (NUM_ROWS / NUM_LAYERS / NUM_GRIDS_X / NUM_GRIDS_Y); 
		int grid_id_x = grid_id / NUM_GRIDS_Y; 
		int grid_id_y = grid_id % NUM_GRIDS_Y; 

		int bank_id_x = bank_id_ / bank_y; 
		int bank_id_y = bank_id_ % bank_y; 

		*row = vault_id_x * (bank_x * NUM_GRIDS_X) + bank_id_x * NUM_GRIDS_X + grid_id_x; 
        *col = vault_id_y * (bank_y * NUM_GRIDS_Y) + bank_id_y * NUM_GRIDS_Y + grid_id_y;
	}
	else if (ARCH_SCHEME == 1)
	{
		// layer # is determined by the index of bank
		// each bank is divided into NUM_GRIDS_X * NUM_GRIDS_Y thermal grids 
		// all the thermal grids within one bank lie on the same layer
		int num_bank_per_layer = NUM_BANKS / NUM_LAYERS; 
		*layer = bank_id_ / num_bank_per_layer; 

		int bank_same_layer = bank_id_ % num_bank_per_layer; 
		int bank_id_x = bank_same_layer / bank_y; 
		int bank_id_y = bank_same_layer % bank_y; 

		int grid_step = NUM_ROWS / (NUM_GRIDS_X * NUM_GRIDS_Y); 
		int grid_id = row_id_ / grid_step; 
		int grid_id_x = grid_id / NUM_GRIDS_Y; 
		int grid_id_y = grid_id % NUM_GRIDS_Y; 

        *row = vault_id_x * (bank_x * NUM_GRIDS_X) + bank_id_x * NUM_GRIDS_X + grid_id_x; 
        *col = vault_id_y * (bank_y * NUM_GRIDS_Y) + bank_id_y * NUM_GRIDS_Y + grid_id_y;
	}
	else if (ARCH_SCHEME == 2)
	{
		int num_bank_per_layer = NUM_BANKS / NUM_LAYERS; 
		*layer = bank_id_ / num_bank_per_layer; 

		*row = vault_id_x;
		*col = vault_id_y;
	}

}

void ThermalCalculator::printP_new(uint64_t cur_cycle){
	genTotalP(true, cur_cycle);
	uint64_t ElapsedCycle = cur_cycle; 
	std::ofstream power_file; 
	power_file.open(avg_power_str.c_str()); 
	power_file << "layer_type,z,x,y,power,vault,bank\n";
	for (int iz = 0; iz < z; iz ++){
		for (int iy = 0; iy < y; iy ++){
			for (int ix = 0; ix < x; ix ++){
				power_file << "MEM," << iz << "," << ix << "," << iy << "," << accu_Pmap_wLogic[ix][iy][iz] / (double) ElapsedCycle / CPU_CLK_PERIOD * tCK << "," << ix/(x/vault_x)*vault_y + iy/(y/vault_y) << "," << (ix%(x/vault_x))/NUM_GRIDS_X*bank_y + (iy%(y/vault_y))/NUM_GRIDS_Y  <<std::endl;
			}
		}
	}
	for (int iy = 0; iy < y; iy ++){
		for (int ix = 0; ix < x; ix ++){
			power_file << "LOGIC," << z << "," << ix << "," << iy << "," << accu_Pmap_wLogic[ix][iy][z] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK << ",-1,-1," << std::endl;
		}
	}
	for (int iy = 0; iy < y; iy ++){
		for (int ix = 0; ix < x; ix ++){
			power_file << "CPU," << z+1 << "," << ix << "," << iy << "," << accu_Pmap_wLogic[ix][iy][z+1] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK << ",-1,-1," << std::endl;
		}
	}
	power_file.close();
}


void ThermalCalculator::printTtrans(unsigned S_id)
{
	// extract the temperature to a 3D array 
	int numP, dimX, dimZ; 
	//double T0 = 273;
	dimX = x; dimZ = y; 
	if (withLogic){
		numP = z + 2;
	}
	else{
		numP = z; 
	}

	double maxTT; maxTT = get_maxT(T_trans, dimX*dimZ*(numP*3+1));
    std::cout << "Now maxT = " << maxTT - T0 << std::endl;

	vector<vector<vector<double> > > T;
	vector<int> layerP; 
	T = vector<vector<vector<double> > > (dimX, vector<vector<double> > (dimZ, vector<double> (numP, 0)));
	layerP = vector<int> (numP, 0);

	/*
	for (int l = 0; l < numP; l ++){
		layerP[l] = l * 3;
        for (int i = 0; i < dimX; i ++){
            for (int j = 0; j < dimZ; j++){
                T[i][j][l] = T_trans[dimX*dimZ*(layerP[l]+1) + j*dimX + i];
            }
        }
	}*/

	for (int l = 0; l < numP; l ++){
		layerP[l] = l * 3;
        for (int i = 0; i < dimX; i ++){
            for (int j = 0; j < dimZ; j++){
                T[i][j][l] = T_trans[dimX*dimZ*(layerP[l]+1) + i*dimZ + j];
            }
        }
	}

	/////////////// print out to files ///////////////
	std::ofstream temp_file; 
	temp_file.open(temp_trace_str.c_str(), std::ios_base::app); 
	for (int iz = 0; iz < numP; iz ++){
		for (int iy = 0; iy < dimZ; iy ++){
			for (int ix = 0; ix < dimX; ix ++){
				temp_file << S_id << "," << iz << "," << ix << "," << iy << "," << T[ix][iy][iz] - T0 << std::endl;
			}
		}
	}
	temp_file.close(); 

}

void ThermalCalculator::printSamplePower2(uint64_t cur_cycle, unsigned S_id){
	uint64_t ElapsedCycle = cur_cycle; 
	std::ofstream power_file; 
	power_file.open(power_trace_str.c_str(), std::ios_base::app);
	if (withLogic){
		for (int iz = 0; iz < z+2; iz ++){
			for (int iy = 0; iy < y; iy ++){
				for (int ix = 0; ix < x; ix ++){
					power_file << S_id << "," << iz << "," << ix << "," << iy << "," << cur_Pmap_wLogic[ix][iy][iz] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK << std::endl;
				}
			}
		}
		power_file.close(); 
	}
	else{
		for (int iz = 0; iz < z; iz ++){
			for (int iy = 0; iy < y; iy ++){
				for (int ix = 0; ix < x; ix ++){
					power_file << S_id << "," << iz << "," << ix << "," << iy << "," << cur_Pmap_wLogic[ix][iy][iz] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK << std::endl;
				}
			}
		}
		power_file.close();
	}
}


void ThermalCalculator::printT()
{
	// print out the temperature profile calcualted using the accumulated power 
	std::ofstream temp_file; 
	temp_file.open(final_temp_str.c_str()); 
	temp_file << "layer,x,y,temperature\n"; 

	int numlayer; 
	if (withLogic)
		numlayer = z + 2; 
	else
		numlayer = z; 

	for (int iz = 0; iz < numlayer; iz ++){
		for (int iy = 0; iy < y; iy ++){
			for (int ix = 0; ix < x; ix ++){
				temp_file << iz << "," << ix << "," << iy << "," << T_final[ix][iy][iz] << std::endl;
			}
		}
	}
	temp_file.close();
}


void ThermalCalculator::printVaultUsage()
{
	std::cout << "Single cell usage:\n";
	for (int i = 0; i < NUM_VAULTS; i++)
		std::cout << vault_usage_single[i] << ", "; 
	std::cout << std::endl;

	std::cout << "Multi cell usage:\n";
	for (int i = 0; i < NUM_VAULTS; i++)
		std::cout << vault_usage_multi[i] << ", "; 
	std::cout << std::endl;
}

void ThermalCalculator::genTotalP(bool accuP, uint64_t cur_cycle)
{
	/* accuP = true: calculate for the accumulative power */
	/* accuP = false: calculate for the transient power */
 
	//std::cout << "\ncome in the genTotalP\n";

	double logicE, cellE_ratio, cellE; 
	uint64_t ElapsedCycle = cur_cycle; 
	std::vector<std::vector<std::vector<double> > > origP; 
	std::vector<std::vector<std::vector<double> > > newP = accu_Pmap_wLogic; // just initilize it 

	std::vector<std::vector<double> > new_logicP_map; 
	new_logicP_map = imresize2D(logicP_map, logicP_x, logicP_y, x, y);


	std::ofstream power_file; 
	power_file.open(debug_power_resize_str.c_str()); 
	power_file << "x,y,power\n";
	for (int i = 0; i < x; i ++){
		for (int j = 0; j < y; j ++){
			power_file << i << "," << j << "," << new_logicP_map[i][j] << std::endl;
		}
	}
	power_file.close();


	//std::cout << "finish imresize2D\n";
	double val = 0.0;
	for (int i = 0; i < x; i ++)
		for (int j = 0; j < y; j++)
			val += new_logicP_map[i][j];

	if (accuP)
	{
		// calculate the accumulative P
		origP = accu_Pmap; 
		logicE = totalEnergy * 1.83; 
		// the ratio is from "Data compression for thermal mitigation in HMC"
	}
	else
	{
		origP = cur_Pmap; 
		logicE = sampleEnergy * 1.83;
	}
	cellE = logicE / x / y; 
	//cellE_ratio = logicE / val;

	for (int l = 0; l < z+2; l ++)
	{
		if (l == z+1)
		{
			for (int i = 0; i < x; i ++)
				for (int j = 0; j < y; j ++)
					newP[i][j][l] = new_logicP_map[i][j] * (double) ElapsedCycle  * CPU_CLK_PERIOD / tCK;
		}
		else if (l == z)
		{
			for (int i = 0; i < x; i ++)
				for (int j = 0; j < y; j ++)
					newP[i][j][l] = cellE; 
		}
		else
		{
			for (int i = 0; i < x; i ++)
				for (int j = 0; j < y; j ++)
					newP[i][j][l] = origP[i][j][l];
		}
	}

	if (accuP)
		accu_Pmap_wLogic = newP; 
	else
		cur_Pmap_wLogic = newP; 
}


void ThermalCalculator::calcT(uint64_t cur_cycle)
{
	// withLogic = 1: calculate with logic layer 
	// withLogic = 0: calculate only the DRAM layers 
	double ***powerM; 
	int i, j, l; // indicators 
	int numP, dimX, dimZ; 
	// uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
	uint64_t ElapsedCycle = cur_cycle; 

	dimX = x; dimZ = y; 
	if (withLogic){
		numP = z + 2; 
		genTotalP(true, cur_cycle);
	}
	else{
		numP = z; 
	}

	if ( !(powerM = (double ***)malloc(dimX * sizeof(double **))) ) printf("Malloc fails for powerM[].\n");
    for (i = 0; i < dimX; i++)
    {
        if ( !(powerM[i] = (double **)malloc(dimZ * sizeof(double *))) ) printf("Malloc fails for powerM[%d][].\n", i);
        for (j = 0; j < dimZ; j++)
        {
            if ( !(powerM[i][j] = (double *)malloc(numP * sizeof(double))) ) printf("Malloc fails for powerM[%d][%d][].\n", i,j);
        }
    }

    if (withLogic)
    {
    	for (i = 0; i < dimX; i ++)
    		for (j = 0; j < dimZ; j ++)
    			for (l = 0; l < numP; l ++)
    				powerM[i][j][l] = accu_Pmap_wLogic[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
    }
    else
    {
    	for (i = 0; i < dimX; i ++)
    		for (j = 0; j < dimZ; j ++)
    			for (l = 0; l < numP; l ++)
    				powerM[i][j][l] = accu_Pmap[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
    }


    double sum_power = 0.0; 
    for (i = 0; i < dimX; i ++){
    	for (j = 0; j < dimZ; j ++){
    		for (l = 0; l < numP; l ++){
    			sum_power += powerM[i][j][l];
    		}
    	}
    }
    cout << "sum_power = " << sum_power << endl;


    /*double grid_size = 200e-6; // randomly choose a value
    double W, L; 
    W = grid_size * (double) dimX; 
    L = grid_size * (double) dimZ; */

    T_final = steady_thermal_solver(powerM, ChipX, ChipZ, numP, dimX, dimZ, Midx, MidxSize);

    // print the final temperature profile 
    printT();

}


void ThermalCalculator::save_sampleP(uint64_t cur_cycle, unsigned S_id)
{
	/* In this method, we save the sampling power to the file */
	/* We also calculate the transient temperature at this time point */
	/* The current transient temperature is stored in T_trans */
	/* Also, the current transient temperature is written to the file */


	genTotalP(false, power_epoch); 
	printSamplePower2(power_epoch, S_id); 

	cout << "\ntime = " << float(clock() - t)/CLOCKS_PER_SEC << " [s]\n";  
	cout << "========= solve for Sample " << S_id << "==== Current time is " << power_epoch * (S_id+1) * CPU_CLK_PERIOD * 1e-9 << "[s] ================\n";

	t = clock(); 
	//bool withLogic = true;

	////// calclate the transient temperature ////////
	calc_trans_T();

	printTtrans(S_id);
	////// calculate the transient PDN //////////////
	//TransPDNsolver();

	/////////////// Update the Dynamic Management Information ///////////////
	cout << "num_refresh = " << num_refresh << endl;
	UpdateRefreshCont(); 


}


void ThermalCalculator::calc_trans_T()
{
	double time = power_epoch * CPU_CLK_PERIOD * 1e-9; // [s]
	// withLogic = 1: calculate with logic layer 
	// withLogic = 0: calculate only the DRAM layers 
	double ***powerM; 
	int i, j, l; // indicators 
	int numP, dimX, dimZ; 
	// uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
	uint64_t ElapsedCycle = power_epoch;

	dimX = x; dimZ = y; 
	if (withLogic){
		numP = z + 2; 
	}
	else{
		numP = z; 
	}

	if ( !(powerM = (double ***)malloc(dimX * sizeof(double **))) ) printf("Malloc fails for powerM[].\n");
    for (i = 0; i < dimX; i++)
    {
        if ( !(powerM[i] = (double **)malloc(dimZ * sizeof(double *))) ) printf("Malloc fails for powerM[%d][].\n", i);
        for (j = 0; j < dimZ; j++)
        {
            if ( !(powerM[i][j] = (double *)malloc(numP * sizeof(double))) ) printf("Malloc fails for powerM[%d][%d][].\n", i,j);
        }
    }

    if (withLogic)
    {
    	for (i = 0; i < dimX; i ++)
    		for (j = 0; j < dimZ; j ++)
    			for (l = 0; l < numP; l ++)
    				powerM[i][j][l] = cur_Pmap_wLogic[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
    }
    else
    {
    	for (i = 0; i < dimX; i ++)
    		for (j = 0; j < dimZ; j ++)
    			for (l = 0; l < numP; l ++)
    				powerM[i][j][l] = cur_Pmap[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
    }

    /*double grid_size = 200e-6; // randomly choose a value
    double W, L; 
    W = grid_size * (double) dimX; 
    L = grid_size * (double) dimZ; */

    int time_iter = TimeIter0; 
    while (time / time_iter >= max_tp)
    	time_iter ++; 


    T_trans = transient_thermal_solver(powerM, ChipX, ChipZ, numP, dimX, dimZ, Midx, MidxSize, Cap, CapSize, time, time_iter, T_trans);

    //double maxTT; maxTT = get_maxT(T_trans, dimX*dimZ*(numP*3+1));
    //std::cout << "Now maxT = " << maxTT - T0 << std::endl;

}

void ThermalCalculator::calcMidx()
{
	int dimX, dimZ, numP;

	dimX = x; dimZ = y; 
	if (withLogic){
		numP = z + 2; 
	}
	else{
		numP = z; 
	}

	/*double grid_size = 200e-6; 
	double W, L; 
	W = grid_size * (double) dimX; 
	L = grid_size * (double) dimZ; */

	layerP_ = vector<int> (numP, 0);
	for(int l = 0; l < numP; l++)
        layerP_[l] = l * 3;



	Midx = calculate_Midx_array(ChipX, ChipZ, numP, dimX, dimZ, &MidxSize);
	Cap = calculate_Cap_array(ChipX, ChipZ, numP, dimX, dimZ, &CapSize);
	calculate_time_step();
	T_trans = initialize_Temperature(ChipX, ChipZ, numP, dimX, dimZ);

    /*cout << "MidxSize = " << MidxSize << endl; 

	int iidx; 
  	for (iidx = 0; iidx < MidxSize; iidx ++)
  		cout << Midx[iidx][0] << "\t" << Midx[iidx][1] << "\t" << Midx[iidx][2] << endl;
*/
}

void ThermalCalculator::calculate_time_step()
{
	double dt = 100.0; 
	int layer_dim = x * y; 
	double c, g; 

	for (int l = 0; l < CapSize; l ++)
	{
		c = Cap[l]; 
		for (int i = 0; i < layer_dim; i ++)
		{
			if (Midx[l*layer_dim + i][0] == Midx[l*layer_dim + i][1])
			{
				g = Midx[l*layer_dim + i][2]; 
				if (c/g < dt) 
					dt = c/g; 
			}
		}
	}

	cout << "maximum dt is " << dt << endl;
	max_tp = dt; 

}


void ThermalCalculator::ReadlogicP()
{
	ifstream filein;
	filein.open(logicPFileName.c_str());
	double totalpower; 

	filein >> logicP_x; 
	filein >> logicP_y; 
	filein >> logicP_z; 
	cout << "logicP_z = " << logicP_z << endl;
	logicP_map = vector<vector<double> > (logicP_x, vector<double> (logicP_y, 0));
	for (int j = 0; j < logicP_y; j ++)
		for (int i = 0; i < logicP_x; i ++)
			filein >> logicP_map[i][j];

	totalpower = 0.0; 
	for (int j = 0; j < logicP_y; j ++)
		for (int i = 0; i < logicP_x; i ++)
			totalpower += logicP_map[i][j];

	cout << "totalpower = " << totalpower << endl;

	std::ofstream power_file; 
	power_file.open(debug_power_str.c_str()); 
	power_file << "x,y,power\n";
	for (int i = 0; i < logicP_x; i ++){
		for (int j = 0; j < logicP_y; j ++){
			power_file << i << "," << j << "," << logicP_map[i][j] << std::endl;
		}
	}
	power_file.close();

}


int ThermalCalculator::get_x()
{
	return x; 
}
int ThermalCalculator::get_y()
{
	return y; 
}
int ThermalCalculator::get_z()
{
	return z; 
}
double ThermalCalculator::get_totalE()
{
	return totalEnergy;
}
double ThermalCalculator::get_IOE()
{
	return IOEnergy;
}

/* Methods for the Dynamic Management */
void ThermalCalculator::UpdateRefreshCont()
{
	int vid, bid, rid_s, rid_e; 
	double T_local; 

	for (int ix = 0; ix < x; ix ++){
		for (int iy = 0; iy < y; iy ++){
			for (int iz = 0; iz < z; iz ++){
				T_local = T_trans[x*y*(layerP_[iz]+1) + iy*x + ix] - T0; 
				//cout << "\rT = " << T_local << flush; 
				rev_mapPhysicalLocation(&vid, &bid, &rid_s, &rid_e, iz, ix, iy);
				RefreshCont.UpdateRetT(vid, bid, rid_s, rid_e, T_local);
			}
		}
	}
	//cout << endl;
}

void ThermalCalculator::printRT(uint64_t cur_cycle)
{
	cout << "Print Retention Time Count Down\n";
	//std::ostringstream file_oss; 
	std::ofstream RT_file;
	std::ofstream cur_file; 

	//file_oss << "./power_trace/RT_sample_" << S_id << ".csv";
	//std::string file_name_str = file_oss.str();
	//file_oss.str("");
	//char* file_name = new char[file_name_str.size() + 1]; 
	//std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
	//file_name[file_name_str.size()] = '\0';
	RT_file.open(dump_RT_str.c_str()); 
	//RT_file << NUM_VAULTS << " " << NUM_BANKS << " " << NUM_ROWS << endl; 

	for (int iv = 0; iv < NUM_VAULTS; iv ++){
		for (int ib = 0; ib < NUM_BANKS; ib ++){
			for (int ir = 0; ir < NUM_ROWS; ir ++){
				RT_file << RefreshCont.RetTCountDown[iv][ib][ir]; 
				if (ir < NUM_ROWS - 1)
					RT_file << " ";
			}
			if (ib < NUM_BANKS - 1)
				RT_file << " ";
		}
		RT_file << ";\n";
	}
	RT_file.close();

	cur_file.open(dump_curCyc_str.c_str());
	cur_file << cur_cycle << endl;
	cur_file.close(); 

}