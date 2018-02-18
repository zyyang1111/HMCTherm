#include "RefreshControl.h"
#include <fstream>
#include <cstdlib>

using namespace std; 
using namespace CasHMC; 

extern string RTFileName;
extern string RT_cntdown_str; 
extern int cont_bool; 


RFControl::RFControl():
	Tref(85.0),
	minRefT(5),
	maxRefT(9)
{
	RetT_ref = vector<vector<vector<double> > > (NUM_VAULTS, vector<vector<double> > (NUM_BANKS, vector<double> (NUM_ROWS, 0)));
	RetT = vector<vector<vector<double> > > (NUM_VAULTS, vector<vector<double> > (NUM_BANKS, vector<double> (NUM_ROWS, 0)));
	RetTCountDown_r = vector<vector<vector<int> > > (NUM_VAULTS, vector<vector<int> > (NUM_BANKS, vector<int> (NUM_ROWS, 0)));
	RetTCountDown = vector<vector<vector<int> > > (NUM_VAULTS, vector<vector<int> > (NUM_BANKS, vector<int> (NUM_ROWS, 0)));


	IniRet(); 
	for (int iv = 0; iv < NUM_VAULTS; iv ++)
		for (int ib = 0; ib < NUM_BANKS; ib ++)
			for (int ir = 0; ir < NUM_ROWS; ir ++)
				RetT[iv][ib][ir] = RetT_ref[iv][ib][ir];

	IniRetCountD(); // calculate this first to get RetTCountDown_r 
					// we will over write RetTCountDown if cont_bool == 1

	if (cont_bool){
		LoadRetCountD();
	}
	
}

void RFControl::LoadRetCountD()
{
	string line; 
	ifstream myfile(RT_cntdown_str.c_str()); 

	cout << "Try to load RetTCountDown_file.txt\n";

	if (myfile.is_open())
	{
		for (int i = 0; i < NUM_VAULTS; i++){
			cout << "\rvault " << i << flush;
			getline(myfile, line);
			parse_line_2(line, i); 
		}
		cout << endl;
		myfile.close();
	}
	else
		cout << "Cannot open RetTCountDown_file.txt\n";
}

void RFControl::parse_line_2(string line, int vaultID)
{
	size_t found;
	int ir = 0, ib = 0;  
	while ((found = line.find_first_of(" ")) != string::npos){
		RetTCountDown[vaultID][ib][ir] = atof(line.substr(0, found).c_str()); 
		line.erase(0, found+1); 
		ir ++; 
		if (ir == NUM_ROWS){
			ir = 0;
			ib ++; 
		}
	}
	found = line.find_first_of(";"); 
	RetTCountDown[vaultID][ib][ir] = atof(line.substr(0, found).c_str()); 
}


void RFControl::IniDim(int x_, int y_, int z_)
{
	x = x_; y = y_; z = z_;
}

void RFControl::IniRet()
{
	string line; 
	ifstream myfile(RTFileName.c_str()); 

	cout << "Try to parse DRAM_RT.txt\n";

	if (myfile.is_open())
	{
		for (int i = 0; i < NUM_VAULTS; i++){
			cout << "\rvault " << i << flush;
			getline(myfile, line);
			parse_line(line, i); 
		}
		cout << endl;
		myfile.close();
	}
	else
		cout << "Cannot open DRAM_RT.txt\n";
}


void RFControl::parse_line(string line, int vaultID)
{
	size_t found;
	int ir = 0, ib = 0;  
	while ((found = line.find_first_of(" ")) != string::npos){
		RetT_ref[vaultID][ib][ir] = atof(line.substr(0, found).c_str()); 
		line.erase(0, found+1); 
		ir ++; 
		if (ir == NUM_ROWS){
			ir = 0;
			ib ++; 
		}
	}
	found = line.find_first_of(";"); 
	RetT_ref[vaultID][ib][ir] = atof(line.substr(0, found).c_str()); 
}

void RFControl::UpdateRetT(int vault, int bank, int row_s, int row_e, int T)
{
	int p; 

	for (int ir = row_s; ir < row_e; ir ++)
	{
		if (RetT_ref[vault][bank][ir] == 0){
			RetT[vault][bank][ir] = 0;
			RetTCountDown_r[vault][bank][ir];
		}
		else{
			RetT[vault][bank][ir] = exp((Tref + T0)/(T + T0) * log(RetT_ref[vault][bank][ir]));
			//cout << "RetT_ref = " << RetT_ref[vault][bank][ir] << "; RetT = " << RetT[vault][bank][ir] << endl;
			// update RetTCountDown_r
			p = log(RetT[vault][bank][ir]) / log(2); 
			//cout << "p = " << p << endl; 
			RetTCountDown_r[vault][bank][ir] = min(exp(log(2) * (p-minRefT)) - 1, (double) maxRefT);
			//cout << "count down = " << RetTCountDown_r[vault][bank][ir] << endl;
		}
	}
}

void RFControl::IniRetCountD()
{
	int p; 
	for (int iv = 0; iv < NUM_VAULTS; iv ++){
		for (int ib = 0; ib < NUM_BANKS; ib ++){
			for (int ir = 0; ir < NUM_ROWS; ir ++){
				if (RetT[iv][ib][ir] ==0){
					RetTCountDown_r[iv][ib][ir] = 0;
				}
				else{
					p = log(RetT[iv][ib][ir]) / log(2); 
					RetTCountDown_r[iv][ib][ir] = min(exp(log(2) * (p-minRefT)) - 1, (double) maxRefT); 
				}
				RetTCountDown[iv][ib][ir] = RetTCountDown_r[iv][ib][ir];
			}
		}
	}
} 

bool RFControl::UpdateCountD(int vault, int bank, int row)
{
	//if (row % 500 == 0)
	//	cout << "RTCD[" << vault << "][" << bank << "][" << row << "] = " << RetTCountDown[vault][bank][row] << endl;
	bool refresh = false; 
	if (!RetTCountDown[vault][bank][row]){
		refresh = true; 
		RetTCountDown[vault][bank][row] = RetTCountDown_r[vault][bank][row]; 
	}
	else
		RetTCountDown[vault][bank][row] --;

	return refresh; 
}

void RFControl::ResetCountD(int vault, int bank, int row)
{
	RetTCountDown[vault][bank][row] = RetTCountDown_r[vault][bank][row]; 
}