#include "Thermal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> //stringstream
#include <stdlib.h> // getenv()


using namespace std; 
using namespace CasHMC; 

extern "C" double ***steady_thermal_solver(double ***powerM, double W, double Lc, int numP, int dimX, int dimZ, double **Midx, int count);
extern "C" double *transient_thermal_solver(double ***powerM, double W, double L, int numP, int dimX, int dimZ, double **Midx, int MidxSize, double *Cap, int CapSize, double time, int iter, double *T_trans);
extern "C" double **calculate_Midx_array(double W, double Lc, int numP, int dimX, int dimZ, int* MidxSize);
extern "C" double *calculate_Cap_array(double W, double Lc, int numP, int dimX, int dimZ, int* CapSize); 
extern "C" double *initialize_Temperature(double W, double Lc, int numP, int dimX, int dimZ); 

ThermalCalculator::ThermalCalculator(bool withLogic_):
	totalEnergy(0.0),
	sampleEnergy(0.0),
	pe_crit(false),
	power_epoch(1000),
	sample_id(0),
	withLogic(withLogic_)
	{
		std::cout << "enter the assignment method\n";
		std::cout << "ARCH_SCHEME = " << ARCH_SCHEME << std::endl;

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
			vault_x = square_array(NUM_VAULTS); 
			bank_y = square_array(num_bank_per_layer); 

			vault_y = NUM_VAULTS / vault_x; 
			bank_x = num_bank_per_layer / bank_y; 
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

	}

ThermalCalculator::~ThermalCalculator()
{
	std::cout << "delete ThermalCalculator\n";
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

void ThermalCalculator::addPower(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, bool single_bank, uint64_t cur_cycle)
{

	//std::cout << "energy = " << energy_t_ << std::endl; 
	//std::cout << "(vault, bank, row, col) = " << "( " << vault_id_ << ", " << bank_id_ << ", " << row_id_ << ", " << col_id_ << " )" << std::endl;
	//std::cout << "single_bank is " << single_bank << std::endl;
    
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

void ThermalCalculator::printP(uint64_t cur_cycle)
{
	std::ostringstream file_oss; 

	//uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
	uint64_t ElapsedCycle = cur_cycle; 

	std::ofstream power_file; 
	for (int iz = 0; iz < z; iz ++)
	{
		file_oss << "power_mem_layer" << iz << ".csv";
		std::string file_name_str = file_oss.str();
		file_oss.str("");
		char* file_name = new char[file_name_str.size() + 1]; 
		std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
		file_name[file_name_str.size()] = '\0';
		power_file.open(file_name); 
		for (int iy = 0; iy < y; iy ++)
		{
			for (int ix = 0; ix < x; ix ++)
			{
				power_file << accu_Pmap[ix][iy][iz] / (double) ElapsedCycle / CPU_CLK_PERIOD * tCK; 
				//power_file << accu_Pmap[ix][iy][iz]; 
				if (ix < x-1)
					power_file << ",";
			}
			power_file << endl; 
		} 
		power_file.close();
	}
}

void ThermalCalculator::printTtrans(unsigned S_id)
{
	// extract the temperature to a 3D array 
	int numP, dimX, dimZ; 
	//double T0 = 273;
	dimX = x; dimZ = y; 
	if (withLogic){
		numP = z + 1;
	}
	else{
		numP = z; 
	}

	vector<vector<vector<double> > > T;
	vector<int> layerP; 
	T = vector<vector<vector<double> > > (dimX, vector<vector<double> > (dimZ, vector<double> (numP, 0)));
	layerP = vector<int> (numP, 0);

	for (int l = 0; l < numP; l ++){
		layerP[l] = l * 3;
        for (int i = 0; i < dimX; i ++){
            for (int j = 0; j < dimZ; j++){
                T[i][j][l] = T_trans[dimX*dimZ*(layerP[l]+1) + j*dimX + i];
            }
        }
	}

	/////////////// print out to files ///////////////

	std::ostringstream file_oss; 
	std::ofstream temp_file;

	file_oss << "./temperature_trace/temp_sample_" << S_id << ".csv"; 
	std::string file_name_str = file_oss.str();
	file_oss.str("");
	char* file_name = new char[file_name_str.size() + 1]; 
	std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
	file_name[file_name_str.size()] = '\0';
	temp_file.open(file_name); 
	for (int iz = 0; iz < numP; iz ++){
		for (int iy = 0; iy < dimZ; iy ++){
			for (int ix = 0; ix < dimX; ix ++){
				temp_file << T[ix][iy][iz] - T0; 
				if (ix < dimX-1)
					temp_file << ","; 
			}
			if (iy < dimZ-1)
				temp_file << ",";
		}
		temp_file << endl; 
	}
	temp_file.close();

}


void ThermalCalculator::printSamplePower(uint64_t cur_cycle, unsigned S_id)
{
	std::ostringstream file_oss; 
    //uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
	uint64_t ElapsedCycle = cur_cycle; 
	std::ofstream power_file;

	file_oss << "./power_trace/power_sample_" << S_id << ".csv"; 
	std::string file_name_str = file_oss.str();
	file_oss.str("");
	char* file_name = new char[file_name_str.size() + 1]; 
	std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
	file_name[file_name_str.size()] = '\0';
	power_file.open(file_name); 

	if (withLogic){
		for (int iz = 0; iz < z+2; iz ++){
			for (int iy = 0; iy < y; iy ++){
				for (int ix = 0; ix < x; ix ++){
					power_file << cur_Pmap_wLogic[ix][iy][iz] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
					if (ix < x-1)
						power_file << ","; 
				}
				if (iy < y-1)
					power_file << ",";
			}
			power_file << endl; 
		}	
		power_file.close();
	}
	else{
		for (int iz = 0; iz < z; iz ++){
			for (int iy = 0; iy < y; iy ++){
				for (int ix = 0; ix < x; ix ++){
					power_file << cur_Pmap[ix][iy][iz] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
					if (ix < x-1)
						power_file << ","; 
				}
				if (iy < y-1)
					power_file << ",";
			}
			power_file << endl; 
		}	
		power_file.close();
	}

}


void ThermalCalculator::printT()
{
	// print out the temperature profile calcualted using the accumulated power 
	std::ostringstream file_oss; 

	std::ofstream temp_file; 

	int numlayer; 
	if (withLogic)
		numlayer = z + 2; 
	else
		numlayer = z; 

	for (int iz = 0; iz < numlayer; iz ++)
	{
		file_oss << "temperature_layer" << iz << ".csv";
		std::string file_name_str = file_oss.str();
		file_oss.str("");
		char* file_name = new char[file_name_str.size() + 1]; 
		std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
		file_name[file_name_str.size()] = '\0';
		temp_file.open(file_name); 
		for (int iy = 0; iy < y; iy ++)
		{
			for (int ix = 0; ix < x; ix ++)
			{
				temp_file << T_final[ix][iy][iz]; 
				if (ix < x - 1)
					temp_file << ",";
			}
			temp_file << endl; 
		}
		temp_file.close();
	}
}


void ThermalCalculator::print_logicP(uint64_t cur_cycle)
{
    genTotalP(true, cur_cycle);

    // uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
    uint64_t ElapsedCycle = cur_cycle;

	std::ofstream power_file; 
	power_file.open("power_processor_layer.csv");
	for (int iy = 0; iy < y; iy ++)
	{
		for (int ix = 0; ix < x; ix ++)
		{
			power_file << accu_Pmap_wLogic[ix][iy][z+1] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
			if (ix < x - 1)
				power_file << ",";
		}
		power_file << endl; 
	}
	power_file.close();


	power_file.open("power_logic_layer.csv");
	for (int iy = 0; iy < y; iy ++)
	{
		for (int ix = 0; ix < x; ix ++)
		{
			power_file << accu_Pmap_wLogic[ix][iy][z] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
			if (ix < x - 1)
				power_file << ",";
		}
		power_file << endl; 
	}
	power_file.close();
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

	std::cout << "come in the genTotalP\n";

	double logicE, cellE_ratio, cellE; 
	uint64_t ElapsedCycle = cur_cycle; 
	std::vector<std::vector<std::vector<double> > > origP; 
	std::vector<std::vector<std::vector<double> > > newP = accu_Pmap_wLogic; // just initilize it 

	std::vector<std::vector<double> > new_logicP_map; 
	new_logicP_map = imresize2D(logicP_map, logicP_x, logicP_y, x, y);

	std::cout << "finish imresize2D\n";
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
	printSamplePower(power_epoch, S_id); 

	cout << "========= solve for Sample " << S_id << "==== Current time is " << power_epoch * (S_id+1) * CPU_CLK_PERIOD * 1e-9 << "[s] ================\n";

	//bool withLogic = true;

	////// calclate the transient temperature ////////
	calc_trans_T();
	printTtrans(S_id);
	////// calculate the transient PDN //////////////
	//TransPDNsolver();


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
	ifstream filein("logicP.in");
	double totalpower; 

	filein >> logicP_x; 
	filein >> logicP_y; 
	logicP_map = vector<vector<double> > (logicP_x, vector<double> (logicP_y, 0));
	for (int j = 0; j < logicP_y; j ++)
		for (int i = 0; i < logicP_x; i ++)
			filein >> logicP_map[i][j];

	totalpower = 0.0; 
	for (int j = 0; j < logicP_y; j ++)
		for (int i = 0; i < logicP_x; i ++)
			totalpower += logicP_map[i][j];

	cout << "totalpower = " << totalpower << endl;

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