#include "Thermal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> //stringstream
#include <stdlib.h> // getenv()
#include <cmath>


using namespace std; 
using namespace CasHMC; 

extern double CPU_CLK_PERIOD;

void ThermalCalculator::IniTransPDN()
{
	V_onC = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	V_offC = vector<double> (2, 0);

	I_onC_x = vector<vector<vector<double> > > (PDN_x-1, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	I_onC_y = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y-1, vector<double> (z+1, 0))); 
	I_onC_z = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	I_offC = vector<double> (2, 0); 

	///// initialize variables ///////
	int i, j, l; 
	// voltage and I_onC_z
	for (i = 0; i < PDN_x; i ++){
		for (j = 0; j < PDN_y; j ++){
			for (l = 0; l < z+1; l ++){
				V_onC[i][j][l] = Vdd; 
				I_onC_z[i][j][l] = 0.0;
			}
		}
	}
	V_offC[0] = V_offC[1] = Vdd; 
	// current 
	for (i = 0; i < PDN_x-1; i ++)
		for (j = 0; j < PDN_y; j ++)
			for (l = 0; l < z+1; l ++)
				I_onC_x[i][j][l] = 0.0;
	for (i = 0; i < PDN_x; i ++)
		for (j = 0; j < PDN_y-1; j ++)
			for (l = 0; l < z+1; l ++)
				I_onC_y[i][j][l] = 0.0;
}


void ThermalCalculator::TransPDNsolver()
{
	////// define dX/dt /////////////
	vector<vector<vector<double> > > dV_onC; 
	vector<vector<vector<double> > > dI_onC_x; 
	vector<vector<vector<double> > > dI_onC_y;
	vector<vector<vector<double> > > dI_onC_z;
	vector<double> dV_offC; 
	vector<double> dI_offC;

	////// assign space to variables //
	dV_onC = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	dV_offC = vector<double> (2, 0); 
	dI_onC_x = vector<vector<vector<double> > > (PDN_x-1, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	dI_onC_y = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y-1, vector<double> (z+1, 0))); 
	dI_onC_z = vector<vector<vector<double> > > (PDN_x, vector<vector<double> > (PDN_y, vector<double> (z+1, 0))); 
	dI_offC = vector<double> (2, 0); 

	///// initialize variables ///////
	int i, j, l; 
	int count = 0;

	/////// generate C4map and TSVmap from PDNC4map and PDNTSVmap //
	vector<vector<int> > C4map, TSVmap; 
	C4map = vector<vector<int> > (PDN_x, vector<int> (PDN_y, 0)); 
	TSVmap = vector<vector<int> > (PDN_x, vector<int> (PDN_y, 0));
    for (i = 0; i < PDN_x; i ++){
        for (j = 0; j < PDN_y; j ++){
            TSVmap[i][j] = PDNTSVmap[count]; 
            count ++;
        }
    }
    count = 0; 
    for (i = 0; i < PDN_x; i ++){
        for (j = 0; j < PDN_y; j ++){
            C4map[i][j] = PDNC4map[count]; 
            count ++;
        }
    }
	/////// generate the resized current power map ////////////////
	genTotalP(false, power_epoch);
	vector<vector<vector<double> > > powerM; 
	powerM = imresize(cur_Pmap_wLogic, PDN_x, PDN_y, z+1);
	uint64_t ElapsedCycle = power_epoch; 
	for (i = 0; i < PDN_x; i ++)
    	for (j = 0; j < PDN_y; j ++)
    		for (l = 0; l < z+1; l ++)
    			powerM[i][j][l] = powerM[i][j][l] / (double) ElapsedCycle  / CPU_CLK_PERIOD * tCK; 

    ///////// initialize capacitance values ///////////////////
    double C_SPC; 
    C_SPC = C_DEN * ChipX * ChipZ / PDN_x / PDN_y; 


	//////// iteratively update the variable ///////
    double time = power_epoch * CPU_CLK_PERIOD * 1e-9; // [s]
    double time_step = 1e-16; // [s]
    long num_step; 
    num_step = (long) (time / time_step) + 1;
    time_step = time / (double) num_step;
    double error_tol = 0.00001;
    double accu_time = 0.0;
   	//int num_step = 200; 
   	//double time_step = time / (double) num_step;

    cout << "time of this period is " << time << endl;
   	cout << "PDN solver time step = " << time_step * 1e6 << "[us]; num_step = " << num_step << endl; 

	double val, cap_onC; 
	double minV, minV_pre; 
	minV_pre = Vdd; 

	for (int k = 0; k < num_step; k ++){
		// calculate dV/dt 
		for (l = 0; l < z+1; l ++){
			for (i = 0; i < PDN_x; i ++){
				for (j = 0; j < PDN_y; j ++){
					val = 0.0; cap_onC = C_BLK + C_SPC; 
					if (i-1 >= 0)
						val += I_onC_x[i-1][j][l]; 
					if (i < PDN_x-1)
						val -= I_onC_x[i][j][l];
					if (j-1 >= 0)
						val += I_onC_y[i][j-1][l];
					if (j < PDN_y-1)
						val -= I_onC_y[i][j][l]; 
					if (l == 0){
						cap_onC += C_BUMP / 2; 
						if (C4map[i][j])
							val += I_onC_z[i][j][l];
					}
					else{
						cap_onC += C_PDNTSV / 2; 
						if (TSVmap[i][j])
							val += I_onC_z[i][j][l];
					}
					if (l < z && TSVmap[i][j]){
						cap_onC += C_PDNTSV / 2; 
						val -= I_onC_z[i][j][l];
					}

					dV_onC[i][j][l] = (val - powerM[i][j][l]/Vdd) / cap_onC; 
					//cout << "val = " << val << "; P = " << powerM[i][j][l]/Vdd << "; cap_onC = " << cap_onC << "; dV_onC = " << dV_onC[i][j][l] << endl;
				}
			}
		}
		dV_offC[0] = (I_offC[0] - I_offC[1]) / C_PCB; 
		val = 0;
		for (i = 0; i < PDN_x; i ++){
			for (j = 0; j < PDN_y; j ++){
				if (C4map[i][j])
					val += I_onC_z[i][j][0]; 
			}
		}
		dV_offC[1] = (I_offC[1] - val) / C_PKG;

		// calculate dI/dt --- x
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x-1; i ++)
				for (j = 0; j < PDN_y; j ++)
					dI_onC_x[i][j][l] = (V_onC[i][j][l] - V_onC[i+1][j][l] - I_onC_x[i][j][l]*R_GRID) / L_GRID; 
		// calculate dI/dt --- y
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x; i ++)
				for (j = 0; j < PDN_y-1; j ++)
					dI_onC_y[i][j][l] = (V_onC[i][j][l] - V_onC[i][j+1][l] - I_onC_y[i][j][l]*R_GRID) / L_GRID; 
		// calculate dI/dt --- z
		for (l = 0; l < z+1; l ++){
			for (i = 0; i < PDN_x; i ++){
				for (j = 0; j < PDN_y; j ++){
					if (l == 0)
						if (C4map[i][j])
							dI_onC_z[i][j][l] = (V_offC[1] - V_onC[i][j][l] - I_onC_z[i][j][l]*R_BUMP) / L_BUMP; 
					else
						if (TSVmap[i][j])
							dI_onC_z[i][j][l] = (V_onC[i][j][l-1] - V_onC[i][j][l] - I_onC_z[i][j][l]*R_PDNTSV) / L_PDNTSV; 
				}
			}
		}
		// calculate dI/dt --- off chip 
		dI_offC[0] = (Vdd - V_offC[0] - I_offC[0]*R_PCB) / L_PCB; 
		dI_offC[1] = (V_offC[0] - V_offC[1] - I_offC[1]*R_PKG) / L_PKG; 

		////////////// update voltage and current ////////////
		// update V
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x; i ++)
				for (j = 0; j < PDN_y; j ++)
					V_onC[i][j][l] += dV_onC[i][j][l] * time_step; 
		V_offC[0] += dV_offC[0] * time_step; 
		V_offC[1] += dV_offC[1] * time_step; 

		// update I --- on-chip x
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x-1; i ++)
				for (j = 0; j < PDN_y; j ++)
					I_onC_x[i][j][l] += dI_onC_x[i][j][l] * time_step; 
		// update I --- on-chip y
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x; i ++)
				for (j = 0; j < PDN_y-1; j ++)
					I_onC_y[i][j][l] += dI_onC_y[i][j][l] * time_step; 	
		// update I --- on-chip z
		for (l = 0; l < z+1; l ++)
			for (i = 0; i < PDN_x; i ++)
				for (j = 0; j < PDN_y; j ++)
					I_onC_z[i][j][l] += dI_onC_z[i][j][l] * time_step; 	
		// update I --- off-chip
		I_offC[0] += dI_offC[0] * time_step;
		I_offC[1] += dI_offC[1] * time_step;

		accu_time += time_step;
		if (accu_time >= time)
			break;

		if (!(k%5000)){
			minV = Vdd;
			for (l = 0; l < z+1; l ++){
				for (i = 0; i < PDN_x; i ++){
					for (j = 0; j < PDN_y; j ++){
						if (V_onC[i][j][l] < minV){
							minV = V_onC[i][j][l];
						}
					}
				}
			}
			if (V_offC[0] < minV) minV = V_offC[0];
			if (V_offC[1] < minV) minV = V_offC[1]; 

			//accu_time += time_step * 5000;
			if (k > 5000){
				//if (minV_pre - minV > 0)
					time_step = time_step / (abs(minV_pre - minV) / error_tol);  
			}

			printf("%d: (%.6f [ns]) minV = %.10f   DminV = %.15f; new time_step = %.6f [ps]\n", k, accu_time * 1e9, minV, minV_pre - minV, time_step * 1e12);

			minV_pre = minV;
		}
	} 

	/*cout << endl; 
	//double minV = Vdd; 
	minV = Vdd;
	for (l = 0; l < z+1; l ++){
		for (i = 0; i < PDN_x; i ++){
			for (j = 0; j < PDN_y; j ++){
				if (V_onC[i][j][l] < minV){
					minV = V_onC[i][j][l];
				}
			}
		}
	}
	if (V_offC[0] < minV) minV = V_offC[0];
	if (V_offC[1] < minV) minV = V_offC[1]; 

	//cout << "minV = " << minV << endl;
	printf("minV = %.10f\n", minV);
*/
/*
	double maxV = 0; 
	for (l = 0; l < z+1; l ++){
		for (i = 0; i < PDN_x; i ++){
			for (j = 0; j < PDN_y; j ++){
				if (V_onC[i][j][l] > maxV){
					maxV = V_onC[i][j][l];
				}
			}
		}
	}
	if (V_offC[0] > maxV) maxV = V_offC[0];
	if (V_offC[1] > maxV) maxV = V_offC[1]; 

	//cout << "minV = " << minV << endl;
	printf("maxV = %.10f\n", maxV);
*/

}