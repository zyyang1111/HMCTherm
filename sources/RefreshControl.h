#ifndef REFRESHCONTROL_H
#define REFRESHCONTROL_H

#include "DRAMConfig.h"
#include "SimConfig.h"
#include "ThermalConfig.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

using namespace std; 

namespace CasHMC
{
	class RFControl
	{
		double Tref; // reference temperature in Celsius
		int minRefT; // minimum Retention Time in 2's power
		int maxRefT; // maximum Retention Time in 2's power
	public:
		vector<vector<vector<double> > > RetT_ref; // [vault, bank, row]
		vector<vector<vector<double> > > RetT; // [vault, bank, row]
		vector<vector<vector<int> > > RetTCountDown_r; // used for reset the count down value
		vector<vector<vector<int> > > RetTCountDown; // used to perform count down


		RFControl();
		~RFControl();
		void UpdateRetT(int vault, int bank, int row_s, int row_e, int T);
		void IniRet();
		void IniRetCountD(); 
		void parse_line(string line, int vaultID);
		bool UpdateCountD(int vault, int bank, int row); 
		void ResetCountD(int vault, int bank, int row);
		// called by other routines with the parsed physical address of the memory 
		// return "true" is the corresponding countdown = 0; return "false" otherwise
		// if countdown == 0, reset countdown value
		// if countdown != 0, countdown --
		void LoadRetCountD();
		void parse_line_2(string line, int vaultID);

	};

}




#endif