/************************************************************************
*  Copyright (c) 2017, 		Zhiyuan Yang
*							Ankur Srivastava
*							Bing Shi 
*							University of Maryland
*							zyyang [at] umd [dot] edu
*  All rights reserved. 
*  
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions are met:
*  
*     * Redistributions of source code must retain the above copyright notice,
*        this list of conditions and the following disclaimer.
*  
*     * Redistributions in binary form must reproduce the above copyright notice,
*        this list of conditions and the following disclaimer in the documentation
*        and/or other materials provided with the distribution.
*  
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
*  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
*  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*********************************************************************************/

#ifndef THERMAL_H
#define THERMAL_H 

#include "DRAMConfig.h"
#include "SimConfig.h"
#include "ThermalConfig.h"
#include <vector>
#include <iostream>
#include <ctime>

namespace CasHMC
{

class ThermalCalculator
{
	bool withLogic; // 1 -- calculated with the logic layer 
					// 0 -- otherwise

	double totalEnergy; 
	double sampleEnergy; 
	double IOEnergy;
	int x, y, z;
	int vault_x, vault_y, bank_x, bank_y; // ease for calculation
	std::vector<std::vector<std::vector<double> > > accu_Pmap;
	std::vector<std::vector<std::vector<double> > > accu_Pmap_wLogic;
	std::vector<std::vector<std::vector<double> > > cur_Pmap;
	std::vector<std::vector<std::vector<double> > > cur_Pmap_wLogic;
	std::vector<int> vault_usage_single;
	std::vector<int> vault_usage_multi; 
	std::vector<int> bank_usage_single; 
	std::vector<int> layerP_;

	int num_refresh;

	// statics
	double totRead_E, totWrite_E, totRef_E, totACT_E, totPre_E, totBack_E; 
	double sapRead_E, sapWrite_E, sapRef_E, sapACT_E, sapPre_E, sapBack_E;


	// temperature variable 
	double **Midx; // Midx storing conductance 
	int MidxSize; 
	double *Cap; // Cap storing the capacitance
	int CapSize; 
	double ***T_final; // final temperature file (-3D)
	double *T_trans; // for storage, I only use a vector to store the transient temperature 

	// logic power parameters 
	//int logicP_x, logicP_y, logicP_z; 
	//std::vector<std::vector<double> > logicP_map; 

	// transient control parameters
	unsigned power_epoch; // power sampling period
	double max_tp; // maximum time step
	unsigned sample_id; // index for the sampling nodes

	// output files string
	std::string power_trace_str; // complete path + file for power trace
	std::string temp_trace_str; // complete path + file for temperature trace
	std::string power_stat_str; // complete path + file for power statics for each trace
	std::string avg_power_str; // complete path + file for average power 
	std::string final_temp_str; // complete path + file for static temperature 
	std::string dump_curCyc_str;
	std::string dump_Ttrans_str; 
	std::string dump_accuP_str; 
	std::string dump_curP_str;
	std::string dump_Pstat_str;

	clock_t t; 

public: 
	ThermalCalculator(bool withLogic_);
	~ThermalCalculator();

    /***********************************************************************/
    /*                   MEHTODS FOR ACQURING THE POWER                    */
    /***********************************************************************/

	void addPower_refresh(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, uint64_t cur_cycle);
	void addPower(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, bool single_bank, uint64_t cur_cycle, int cmd_type); 
	//void addIOPower(double energy_t_, unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, uint64_t cur_cycle);
	int square_array(int total_grids_);
	void mapPhysicalLocation(unsigned vault_id_, unsigned bank_id_, unsigned row_id_, unsigned col_id_, int *layer, int *row, int *col);
    void rev_mapPhysicalLocation(int *vault_id_, int *bank_id_, int *row_s, int *row_e, int layer, int row, int col);
    
    // for determining the size of the power map by the external methods
	int get_x(); 
	int get_y();
	int get_z();
	// for acquring the total power by the external methods
	double get_totalE();
	double get_IOE();
	// print accu_Pmap 
	//void printP(uint64_t cur_cycle); 
	void printP_new(uint64_t cur_cycle); 
	//void print_logicP(uint64_t cur_cycle);
	void printVaultUsage();
	void genTotalP(bool accuP, uint64_t cur_cycle);
	void Dump_PTdata();
	void Reload_PTdata();
	void printRT(uint64_t cur_cycle);

	std::vector<std::vector<std::vector<double> > > imresize(std::vector<std::vector<std::vector<double> > > InImage, int x_new, int y_new, int z_);
    std::vector<std::vector<double> > imresize2D(std::vector<std::vector<double> > InImage, int x_old, int y_old, int x_new, int y_new);

    /***********************************************************************/
    /*                MEHTODS FOR PROCESSING THE TEMPERATURE               */
    /***********************************************************************/

	// calculate temperature
	void calcT(uint64_t cur_cycle);
	void calcMidx();
	void calc_trans_T();

	// print the steady state temperature 
	void printT();

	// process the transient thermal simulation 
	void save_sampleP(uint64_t cur_cycle, unsigned S_id);
	void printSamplePower(uint64_t cur_cycle, unsigned S_id); // S_id -> sampling id of the current feature map
	void printTtrans(unsigned S_id); 
	void calculate_time_step();


    /***********************************************************************/
    /*                   Any Dynamic Management Modules                    */
    /***********************************************************************/
    //RFControl RefreshCont; 
    //void UpdateRefreshCont(); 
    

};

}




#endif // THERMAL_H