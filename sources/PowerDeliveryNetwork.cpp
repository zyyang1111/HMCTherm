#include "Thermal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> //stringstream
#include <stdlib.h> // getenv()


using namespace std; 
using namespace CasHMC; 

extern "C" double ***calculate_steady_PDN(double ***powerM, double **Midx, int count, int Wnum, int Lnum, int layer_num, double svdd);
extern "C" double **calculate_PDN_Midx(int **C4map, int **TSVmap, int Wnum, int Lnum, int layer_num, int *MidxSize);

extern double CPU_CLK_PERIOD;

void ThermalCalculator::calcPDNMidx()
{
    ReadC4map(); 
    ReadTSVmap();
    int i, j; // indicators 
    int numP = z+logicP_z+1; 

    ///////////////////// generate the PDN TSV map ////////////////////////////////////////
    int **TSV_map; 
    if ( !(TSV_map = (int **)malloc(PDN_x * sizeof(int *))) ) printf("Malloc fails for TSV_map[].\n");
    for (i = 0; i < PDN_x; i++)
    {
        if ( !(TSV_map[i] = (int *)malloc(PDN_y * sizeof(int))) ) printf("Malloc fails for TSV_map[%d][].\n", i);
    }
    int count = 0; 
    for (i = 0; i < PDN_x; i ++){
        for (j = 0; j < PDN_y; j ++){
            TSV_map[i][j] = PDNTSVmap[count]; 
            count ++;
        }
    }

    ///////////////////// generate the PDN C4 map ////////////////////////////////////////
    int **C4_map; 
    if ( !(C4_map = (int **)malloc(PDN_x * sizeof(int *))) ) printf("Malloc fails for C4_map[].\n");
    for (i = 0; i < PDN_x; i++)
    {
        if ( !(C4_map[i] = (int *)malloc(PDN_y * sizeof(int))) ) printf("Malloc fails for C4_map[%d][].\n", i);
    }
    count = 0; 
    for (i = 0; i < PDN_x; i ++){
        for (j = 0; j < PDN_y; j ++){
            C4_map[i][j] = PDNC4map[count]; 
            count ++;
        }
    }


    PDNMidx = calculate_PDN_Midx(C4_map, TSV_map, PDN_x, PDN_y, numP, &PDNMidxSize);

}


void ThermalCalculator::calc_steadyPDN(uint64_t cur_cycle) 
{
    genTotalP(true, cur_cycle);

	vector<vector<vector<double> > > ResizedP; 
	ResizedP = imresize(accu_Pmap_wLogic, PDN_x, PDN_y, z+logicP_z+1);
    printResizedP(ResizedP, PDN_x, PDN_y, z+logicP_z+1, cur_cycle);

	/////////////////////// generate the power map /////////////////////////////////////////
	double ***powerM; 
	int i, j, l; // indicators 
	int numP, dimX, dimZ; 
	// uint64_t ElapsedCycle = (cur_cycle % LOG_EPOCH == 0)?(LOG_EPOCH):(cur_cycle % LOG_EPOCH); 
	uint64_t ElapsedCycle = cur_cycle; 

	numP = z+logicP_z+1; dimX = PDN_x; dimZ = PDN_y;
	if ( !(powerM = (double ***)malloc(dimX * sizeof(double **))) ) printf("Malloc fails for powerM[].\n");
    for (i = 0; i < dimX; i++)
    {
        if ( !(powerM[i] = (double **)malloc(dimZ * sizeof(double *))) ) printf("Malloc fails for powerM[%d][].\n", i);
        for (j = 0; j < dimZ; j++)
        {
            if ( !(powerM[i][j] = (double *)malloc(numP * sizeof(double))) ) printf("Malloc fails for powerM[%d][%d][].\n", i,j);
        }
    }

    for (i = 0; i < dimX; i ++)
    	for (j = 0; j < dimZ; j ++)
    		for (l = 0; l < numP; l ++)
    			powerM[i][j][l] = ResizedP[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 



    //////////////////// perform the calculation ////////////////////////////////////////
    cout << "======================== PROCESS THE PDN VOLTAGE ============================\n";

    V_final = calculate_steady_PDN(powerM, PDNMidx, PDNMidxSize, dimX, dimZ, numP, Vdd);

    cout << "get V_final \n";
    //// print out the 3D voltage map /////
    printV();

    cout << "finish print!\n";


}

void ThermalCalculator::ReadTSVmap()
{
	ifstream filein("PDN_TSV.map");

    int px, py; 
    filein >> px; 
    filein >> py;

	PDNTSVmap = vector<int> (PDN_x * PDN_y, 0); 
	for (int i = 0; i < PDN_x * PDN_y; i ++)
		filein >> PDNTSVmap[i]; 
}

void ThermalCalculator::ReadC4map()
{
	ifstream filein("PDN_C4.map");

	filein >> PDN_x; 
	filein >> PDN_y; 
	PDNC4map = vector<int> (PDN_x * PDN_y, 0); 
	for (int i = 0; i < PDN_x * PDN_y; i ++)
		filein >> PDNC4map[i]; 
}

void ThermalCalculator::printV()
{
    // print out the voltage profile calcualted using the accumulated power 
    std::ostringstream file_oss; 

    std::ofstream volt_file; 

    for (int iz = 0; iz < z+1; iz ++)
    {
        file_oss << "voltage_layer" << iz << ".csv";
        std::string file_name_str = file_oss.str();
        file_oss.str("");
        char* file_name = new char[file_name_str.size() + 1]; 
        std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
        file_name[file_name_str.size()] = '\0';
        volt_file.open(file_name); 
        for (int iy = 0; iy < PDN_y; iy ++)
        {
            for (int ix = 0; ix < PDN_x; ix ++)
            {
                volt_file << V_final[ix][iy][iz]; 
                if (ix < PDN_x - 1)
                    volt_file << ",";
            }
            volt_file << endl; 
        }
        volt_file.close();
    }
}


void ThermalCalculator::printResizedP(vector<vector<vector<double> > > P, int dimX, int dimY, int dimZ, uint64_t cur_cycle)
{
    // print out the resized profile
    uint64_t ElapsedCycle = cur_cycle;
    std::ostringstream file_oss; 

    std::ofstream RP_file; 

    for (int iz = 0; iz < dimZ; iz ++)
    {
        file_oss << "resizedP_layer" << iz << ".csv";
        std::string file_name_str = file_oss.str();
        file_oss.str("");
        char* file_name = new char[file_name_str.size() + 1]; 
        std::copy(file_name_str.begin(), file_name_str.end(), file_name); 
        file_name[file_name_str.size()] = '\0';
        RP_file.open(file_name); 
        for (int iy = 0; iy < dimY; iy ++)
        {
            for (int ix = 0; ix < dimX; ix ++)
            {
                RP_file << P[ix][iy][iz] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 
                if (ix < dimX - 1)
                    RP_file << ",";
            }
            RP_file << endl; 
        }
        RP_file.close();
    }
}