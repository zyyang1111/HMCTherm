#include "Thermal.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream> //stringstream
#include <stdlib.h> // getenv()


using namespace std; 
using namespace CasHMC; 


vector<vector<vector<double> > > ThermalCalculator::imresize(vector<vector<vector<double> > > InImage, int x_new, int y_new, int z_)
{
	//////////// resize a 3D power map in 2D ////////////////////////////
	/////// The number of layers of the power map keeps unchanged ///////
	////// the order of the layer is flipped ///////////////////////////

	vector<vector<vector<double> > > OutImage; 
	OutImage = vector<vector<vector<double> > > (x_new, vector<vector<double> > (y_new, vector<double> (z_, 0)));
	// z_ is the 3rd dimension of the input image 
	double gridx_in, gridy_in, gridx_out, gridy_out, gridA_in, gridA_out; 

	gridx_in = ChipX/x; gridy_in = ChipZ/y;
	gridx_out = ChipX/x_new; gridy_out = ChipZ/y_new; 

	gridA_in = gridx_in * gridy_in;
	gridA_out = gridx_out * gridy_out;  

	int l, i, j;
	vector<double> total_val; 
	total_val = vector<double> (z_, 0); 
	////////////// calculate the power density ///////////////
	for (l = 0; l < z_; l ++){
		for (i = 0; i < x; i ++){
			for (j = 0; j < y; j ++){
				total_val[l] += InImage[i][j][l];
				InImage[i][j][l] = InImage[i][j][l] / gridA_in; 
			}
		}
	}


	////////////// get the new power density /////////////////
	///////////// using nearest method to do /////////////////

	//printf("x_new = %d; y_new = %d\n", x_new, y_new);

	double gx, gy; 
	int ind_x, ind_y; 
	for (l = 0; l < z_; l ++){
		for (i = 0; i < x_new; i ++){
			for (j = 0; j < y_new; j ++){
				gx = (i+0.5) * gridx_out; gy = (j+0.5) * gridy_out; 
				ind_x = (int) (gx / gridx_in); ind_y = (int) (gy / gridy_in); 
				if (l == 0){
					//printf("(%d, %d): gx(%.2f), gxin(%.2f), indx(%d), gy(%.2f), gyin(%.2f), indy(%d)\n", i,j,gx, gridx_in, ind_x, gy, gridy_in, ind_y);
				}
				OutImage[i][j][l] = InImage[ind_x][ind_y][z_-l-1] * gridA_out;
			}
		}
	}


	//////////// calculate the power for each output grid /////////
	double new_val; 
	for (l = 0; l < z_; l ++){
		new_val = 0.0;
		for (i = 0; i < x_new; i++){
			for (j = 0; j < y_new; j ++){
				new_val += OutImage[i][j][l]; 
			}
		}

		for (i = 0; i < x_new; i ++){
			for (j = 0; j < y_new; j++){
				OutImage[i][j][l] = OutImage[i][j][l] / new_val * total_val[z_-l-1]; 
			}
		}
	}

	return OutImage; 

} 



vector<vector<double> > ThermalCalculator::imresize2D(vector<vector<double> > InImage, int x_old, int y_old, int x_new, int y_new)
{
	//////////// resize a 3D power map in 2D ////////////////////////////
	/////// The number of layers of the power map keeps unchanged ///////
	////// the order of the layer is flipped ///////////////////////////

	vector<vector<double> > OutImage; 
	OutImage = vector<vector<double> > (x_new, vector<double> (y_new, 0));
	// z_ is the 3rd dimension of the input image 
	double gridx_in, gridy_in, gridx_out, gridy_out, gridA_in, gridA_out; 

	gridx_in = ChipX/x_old; gridy_in = ChipZ/y_old;
	gridx_out = ChipX/x_new; gridy_out = ChipZ/y_new; 

	gridA_in = gridx_in * gridy_in;
	gridA_out = gridx_out * gridy_out;  

	int l, i, j;
	double total_val = 0.0; 
	for (i = 0; i < x_old; i ++){
		for (j = 0; j < y_old; j ++){
			total_val += InImage[i][j];
			InImage[i][j] = InImage[i][j] / gridA_in; 
		}
	}

	cout << "total_val = " << total_val << "\n";

	////////////// get the new power density /////////////////
	///////////// using nearest method to do /////////////////

	//printf("x_new = %d; y_new = %d\n", x_new, y_new);

	double gx, gy; 
	int ind_x, ind_y; 
	for (i = 0; i < x_new; i ++){
		for (j = 0; j < y_new; j ++){
			gx = (i+0.5) * gridx_out; gy = (j+0.5) * gridy_out; 
			ind_x = (int) (gx / gridx_in); ind_y = (int) (gy / gridy_in); 
			OutImage[i][j] = InImage[ind_x][ind_y] * gridA_out;
		}
	}



	//////////// calculate the power for each output grid /////////
	double new_val = 0.0; 

	for (i = 0; i < x_new; i++){
		for (j = 0; j < y_new; j ++){
			new_val += OutImage[i][j]; 
		}
	}

	for (i = 0; i < x_new; i ++){
		for (j = 0; j < y_new; j++){
			OutImage[i][j] = OutImage[i][j] / new_val * total_val; 
		}
	}


	return OutImage; 

} 