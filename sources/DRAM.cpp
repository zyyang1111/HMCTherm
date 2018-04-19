/*********************************************************************************
*  CasHMC v1.2 - 2016.09.27
*  A Cycle-accurate Simulator for Hybrid Memory Cube
*
*  Copyright (c) 2016, Dong-Ik Jeon
*                      Ki-Seok Chung
*                      Hanyang University
*                      estwings57 [at] gmail [dot] com
*  All rights reserved.
*********************************************************************************/

#include "DRAM.h"
#include "VaultController.h"

namespace CasHMC
{
	
DRAM::DRAM(ofstream &debugOut_, ofstream &stateOut_, unsigned id, VaultController *contP, ThermalCalculator *tcPtr):
	SimulatorObject(debugOut_, stateOut_),
	DRAMID(id),
	vaultContP(contP),
	prev_state_cycle(0)
{
	classID << DRAMID;
	header = "        (DR_" + classID.str() + ")";
	
	readData = NULL;
	dataCyclesLeft = 0;

	bankStates.reserve(NUM_BANKS);
	for(int b=0; b<NUM_BANKS; b++) {
		bankStates.push_back(new BankState(b));
	}
	previousBankState = vector<BankStateType>(NUM_BANKS, IDLE);

	//thermalCalPtr = new ThermalCalculator();
	thermalCalPtr = tcPtr;
}

DRAM::~DRAM()
{
	// std::cout << "DRAM: currentClockCycle = " << currentClockCycle << std::endl;
	for(int i=0; i<readReturnDATA.size(); i++) {
		delete readReturnDATA[i];
	}
	readReturnDATA.clear();
	readReturnCountdown.clear();
	
	for(int b=0; b<NUM_BANKS; b++) {
		delete bankStates[b];
	}
	bankStates.clear();
}

//
//Receive command from vault controller
//
void DRAM::receiveCMD(DRAMCommand *recvCMD)
{
    int layer, x, y; // 3D physial locations
    unsigned vault_id, bank_id, row_id, col_id, row_id_refresh; // memory address, readable convenience 
    double energy_t = 0.0;
    int refresh_epoch;

    vault_id = vaultContP->vaultContID; 
    bank_id = recvCMD->bank;
    row_id = recvCMD->row; 
    col_id = recvCMD->column;

    //std::cout << "currentClockCycle = " << currentClockCycle << std::endl;

    //std::cout << "current clock cycle = " << currentClockCycle << ": ( " << vault_id << ", " << bank_id << ", " << row_id << ", " << col_id << " )" << std::endl;  
    //std::cout << currentClockCycle << ": vault = " << vault_id << std::endl;

	switch(recvCMD->commandType) {
		case ACTIVATE:
		    energy_t = (((double)IDD0 * (double)tRC) - (((double)IDD3N * (double)tRAS) + ((double)IDD2N * ((double)tRC - (double)tRAS))));
		    thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 1); 
		    //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer

			bankStates[recvCMD->bank]->currentBankState = ROW_ACTIVE;
			bankStates[recvCMD->bank]->lastCommand = ACTIVATE;
			bankStates[recvCMD->bank]->openRowAddress = recvCMD->row;
			bankStates[recvCMD->bank]->nextActivate = max(currentClockCycle + tRC, bankStates[recvCMD->bank]->nextActivate);
			bankStates[recvCMD->bank]->nextPrecharge = max(currentClockCycle + tRAS, bankStates[recvCMD->bank]->nextPrecharge);
			bankStates[recvCMD->bank]->nextRead = max(currentClockCycle + (tRCD-AL), bankStates[recvCMD->bank]->nextRead);
			bankStates[recvCMD->bank]->nextWrite = max(currentClockCycle + (tRCD-AL), bankStates[recvCMD->bank]->nextWrite);
			for(int b=0; b<NUM_BANKS; b++) {
				if(recvCMD->bank != b) {
					bankStates[b]->nextActivate = max(currentClockCycle + tRRD, bankStates[b]->nextActivate);
				}
			}
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      next READ or WRITE time : "<<bankStates[recvCMD->bank]->nextRead<<" [HMC clk]");
			delete recvCMD;
			break;
		case READ:
            energy_t = ((double)IDD4R - (double)IDD3N) * (double)BL/2.0;
            thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 2); 
            energy_t = IOREAD * 32 * 8; // read IO power, spread it across the layer
            thermalCalPtr->addIOPower(energy_t, vault_id, bank_id, row_id, col_id, currentClockCycle);
            //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer
            

			bankStates[recvCMD->bank]->nextPrecharge = max(currentClockCycle + READ_TO_PRE_DELAY, bankStates[recvCMD->bank]->nextPrecharge);
			bankStates[recvCMD->bank]->lastCommand = READ;
			for(int b=0; b<NUM_BANKS; b++) {
				bankStates[b]->nextRead = max(currentClockCycle + max(tCCD, BL), bankStates[b]->nextRead);
				bankStates[b]->nextWrite = max(currentClockCycle + READ_TO_WRITE_DELAY, bankStates[b]->nextWrite);
			}
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      READ DATA return time : "<<currentClockCycle+RL+BL<<" [HMC clk]"
													<<" / next PRECHARGE time : "<<bankStates[recvCMD->bank]->nextPrecharge<<" [HMC clk]");
			recvCMD->commandType = READ_DATA;
			readReturnDATA.push_back(recvCMD);
			readReturnCountdown.push_back(RL);
			break;
		case READ_P:
            energy_t = ((double)IDD4R - (double)IDD3N) * (double)BL/2.0;
            thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 2); 
            energy_t = IOREAD * 32 * 8; // read IO power  
            thermalCalPtr->addIOPower(energy_t, vault_id, bank_id, row_id, col_id, currentClockCycle);
            //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer

			bankStates[recvCMD->bank]->nextActivate = max(currentClockCycle + READ_AUTOPRE_DELAY, bankStates[recvCMD->bank]->nextActivate);
			bankStates[recvCMD->bank]->lastCommand = READ_P;
			bankStates[recvCMD->bank]->stateChangeCountdown = READ_TO_PRE_DELAY;
			for(int b=0; b<NUM_BANKS; b++) {
				bankStates[b]->nextRead = max(currentClockCycle + max(tCCD, BL), bankStates[b]->nextRead);
				bankStates[b]->nextWrite = max(currentClockCycle + READ_TO_WRITE_DELAY, bankStates[b]->nextWrite);
			}
			if(recvCMD->commandType == READ_P) {
				bankStates[recvCMD->bank]->nextRead = bankStates[recvCMD->bank]->nextActivate;
				bankStates[recvCMD->bank]->nextWrite = bankStates[recvCMD->bank]->nextActivate;
			}	
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      READ DATA return time : "<<currentClockCycle+RL+BL<<" [HMC clk]"
													<<" / next ACTIVATE time : "<<bankStates[recvCMD->bank]->nextActivate<<" [HMC clk]");
			recvCMD->commandType = READ_DATA;
			readReturnDATA.push_back(recvCMD);
			readReturnCountdown.push_back(RL);
			break;
		case WRITE:
            energy_t = ((double)IDD4W - (double)IDD3N) * (double)BL/2.0;
            thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 3); 
            energy_t = IOWRITE * 32 * 8; // write IO power  
            thermalCalPtr->addIOPower(energy_t, vault_id, bank_id, row_id, col_id, currentClockCycle);
            //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer

			bankStates[recvCMD->bank]->nextPrecharge = max(currentClockCycle + WRITE_TO_PRE_DELAY, bankStates[recvCMD->bank]->nextPrecharge);
			bankStates[recvCMD->bank]->lastCommand = WRITE;
			for(int b=0; b<NUM_BANKS; b++) {
				bankStates[b]->nextWrite = max(currentClockCycle + max(BL, tCCD), bankStates[b]->nextWrite);
				bankStates[b]->nextRead = max(currentClockCycle + WRITE_TO_READ_DELAY_B, bankStates[b]->nextRead);
			}
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      WRTIE DATA time : "<<currentClockCycle+WL+BL<<" [HMC clk]"
													<<" / next PRECHARGE time : "<<bankStates[recvCMD->bank]->nextPrecharge<<" [HMC clk]");
			delete recvCMD;
			break;
		case WRITE_P:
            energy_t = ((double)IDD4W - (double)IDD3N) * (double)BL/2.0;
            thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 3); 
            energy_t = IOWRITE * 32 * 8; // write IO power  
            thermalCalPtr->addIOPower(energy_t, vault_id, bank_id, row_id, col_id, currentClockCycle);
            //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer

			bankStates[recvCMD->bank]->nextActivate = max(currentClockCycle + WRITE_AUTOPRE_DELAY, bankStates[recvCMD->bank]->nextActivate);
			bankStates[recvCMD->bank]->lastCommand = WRITE_P;
			bankStates[recvCMD->bank]->stateChangeCountdown = WRITE_TO_PRE_DELAY;
			for(int b=0; b<NUM_BANKS; b++) {
				bankStates[b]->nextWrite = max(currentClockCycle + max(BL, tCCD), bankStates[b]->nextWrite);
				bankStates[b]->nextRead = max(currentClockCycle + WRITE_TO_READ_DELAY_B, bankStates[b]->nextRead);
			}
			bankStates[recvCMD->bank]->nextRead = bankStates[recvCMD->bank]->nextActivate;
			bankStates[recvCMD->bank]->nextWrite = bankStates[recvCMD->bank]->nextActivate;
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      WRTIE DATA time : "<<currentClockCycle+WL+BL<<" [HMC clk]"
													<<" / next ACTIVATE time : "<<bankStates[recvCMD->bank]->nextActivate<<" [HMC clk]");
			delete recvCMD; 
			break;
		case WRITE_DATA:
			delete recvCMD;
			break;
		case PRECHARGE:
            energy_t = ((double)IDD0 - (double)IDD2N) * ((double)tRC - (double)tRAS);
            thermalCalPtr->addPower(energy_t, vault_id, bank_id, row_id, col_id, true, currentClockCycle, 4); 
            //thermalCalPtr->mapPhysicalLocation(vault_id, bank_id, row_id, col_id, &layer, &x, &y); // adding the delay of layer

			bankStates[recvCMD->bank]->currentBankState = PRECHARGING;
			bankStates[recvCMD->bank]->lastCommand = PRECHARGE;
			bankStates[recvCMD->bank]->openRowAddress = 0;
			bankStates[recvCMD->bank]->stateChangeCountdown = tRP;
			bankStates[recvCMD->bank]->nextActivate = max(currentClockCycle + tRP, bankStates[recvCMD->bank]->nextActivate);
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      next ACTIVATE time : "<<bankStates[recvCMD->bank]->nextActivate<<" [HMC clk]");
			delete recvCMD;
			break;	
		case REFRESH:
			//cout << "Hey I am here in REFRESH!\n";
			refresh_epoch = (long)((currentClockCycle-100)/(REFRESH_PERIOD/tCK))% (NUM_ROWS / REFRESH_ROWNUM);
			row_id_refresh = refresh_epoch * REFRESH_ROWNUM; 
			//cout << "\nvault_id = " << vault_id << endl; 
			cout << "\rrow_id_refresh = " << row_id_refresh; 

            energy_t = ((double)IDD5 - (double)IDD3N) * (double)tRFC;

			for(int b=0; b<NUM_BANKS; b++) {
				/*
				bool refresh_flag = false; 
				for (int ir = row_id_refresh; ir < row_id_refresh+REFRESH_ROWNUM; ir ++){
					if (thermalCalPtr->RefreshCont.UpdateCountD(vault_id, b, ir)){
						//cout << "refresh!\n";
						refresh_flag = true; 
						//cout << "Issue refresh\n";
						break;
					}
				}
				if (refresh_flag){
					for (int ir = row_id_refresh; ir < row_id_refresh+REFRESH_ROWNUM; ir ++){
						thermalCalPtr->RefreshCont.ResetCountD(vault_id, b, ir);
					}
					thermalCalPtr->addPower_refresh(energy_t, vault_id, b, row_id_refresh, 0, currentClockCycle); 

					bankStates[b]->nextActivate = currentClockCycle + tRFC;
					bankStates[b]->currentBankState = REFRESHING;
					bankStates[b]->lastCommand = REFRESH;
					bankStates[b]->stateChangeCountdown = tRFC;
				}
				*/


				
				thermalCalPtr->addPower_refresh(energy_t, vault_id, b, row_id_refresh, 0, currentClockCycle); 

				bankStates[b]->nextActivate = currentClockCycle + tRFC;
				bankStates[b]->currentBankState = REFRESHING;
				bankStates[b]->lastCommand = REFRESH;
				bankStates[b]->stateChangeCountdown = tRFC;
	            


				//thermalCalPtr->addPower(energy_t, vault_id, b, 0, 0, false, currentClockCycle); 

				//bankStates[b]->nextActivate = currentClockCycle + tRFC;
				//bankStates[b]->currentBankState = REFRESHING;
				//bankStates[b]->lastCommand = REFRESH;
				//bankStates[b]->stateChangeCountdown = tRFC;
			}
			DEBUG(ALI(18)<<header<<ALI(15)<<*recvCMD<<"      next ACTIVATE time : "<<bankStates[recvCMD->bank]->nextActivate<<" [HMC clk]");
			delete recvCMD;
			break;
		default:
			ERROR(header<<"  (DR) == Error - Wrong command type    popped command : "<<*recvCMD<<"  (CurrentClock : "<<currentClockCycle<<")");
			exit(0);
			break;
	}

    //if (recvCMD->commandType != REFRESH)
     //std::cout << "current clock cycle = " << currentClockCycle << ": ( " << vault_id << ", " << bank_id << ", " << row_id << ", " << col_id << " )" << "Energy = " << energy_t << " " << recvCMD->commandType << std::endl; 
	//std::cout << "receive CMD: Energy_t = " << energy_t << std::endl;
}

//
//Switch DRAM to power down mode
//
bool DRAM::powerDown()
{
	//cout << "vault " << vaultContP->vaultContID << ": Power Down -- " << currentClockCycle << endl;
	unsigned layer, x, y; // 3D physial locations
    unsigned vault_id; // memory address, readable convenience 
    double energy_t; 
    uint64_t elapse_time; 

    vault_id = vaultContP->vaultContID; 

	//check to make sure all banks are idle
	bool allIdle = true;
	for(int b=0; b<NUM_BANKS; b++) {
		if(bankStates[b]->currentBankState != IDLE) {
			allIdle = false;
			break;
		}
	}
	if(allIdle) {
		elapse_time = currentClockCycle - prev_state_cycle; 
		prev_state_cycle = currentClockCycle; 
		//if (vault_id == 1)
			//cout << "Vault 1: Power Down -- time = " << elapse_time << "; cur_cycle = " << currentClockCycle << endl;

        energy_t = energy_t = (double)IDD2P * (double)elapse_time / NUM_GRIDS_X / NUM_GRIDS_Y/ NUM_BANKS; // previous state is powerUp
        //std::cout << "powerDown: Energy_t = " << energy_t << std::endl; 

		for(int b=0; b<NUM_BANKS; b++) {
            thermalCalPtr->addPower(energy_t, vault_id, b, 0, 0, false, currentClockCycle, 0); 

			bankStates[b]->currentBankState = POWERDOWN;
			bankStates[b]->nextPowerUp = currentClockCycle + tCKE;
		}
		return true;
	}
	else {
		return false;
	}
}

//
//Awake DRAM in power down mode
//
void DRAM::powerUp()
{
	//cout << "vault " << vaultContP->vaultContID << ": Power Up -- " << currentClockCycle << endl;
	unsigned layer, x, y; // 3D physial locations
    unsigned vault_id; // memory address, readable convenience 
    double energy_t; 
    uint64_t elapse_time;

    vault_id = vaultContP->vaultContID; 

    elapse_time = currentClockCycle - prev_state_cycle; 
    prev_state_cycle = currentClockCycle; 
    //if (vault_id == 1)
		//cout << "Vault 1: Power Up -- time = " << elapse_time << "; cur_cycle = " << currentClockCycle << endl;

    energy_t = (double)IDD3N * (double)elapse_time / NUM_GRIDS_X / NUM_GRIDS_Y / NUM_BANKS;
    //std::cout << "powerUp: Energy_t = " << energy_t << std::endl; 

	for(int b=0; b<NUM_BANKS; b++) {
        thermalCalPtr->addPower(energy_t, vault_id, b, 0, 0, false, currentClockCycle, 0); 

		bankStates[b]->lastCommand = POWERDOWN_EXIT;
		bankStates[b]->currentBankState = AWAKING;
		bankStates[b]->stateChangeCountdown = tXP;
		bankStates[b]->nextActivate = currentClockCycle + tXP;
	}
	DEBUG(ALI(39)<<header<<"AWAKE DRAM power down mode   nextActivate time : "<<bankStates[0]->nextActivate);
}

//
//Updates the state of DRAM
//
void DRAM::Update()
{
	//update bank states
	for(int b=0; b<NUM_BANKS; b++) {
		if(bankStates[b]->stateChangeCountdown>0) {
			//decrement counters
			bankStates[b]->stateChangeCountdown--;

			//if counter has reached 0, change state
			if(bankStates[b]->stateChangeCountdown == 0) {
				switch(bankStates[b]->lastCommand) {
					//only these commands have an implicit state change
					case WRITE_P:
					case READ_P:
						bankStates[b]->currentBankState = PRECHARGING;
						bankStates[b]->lastCommand = PRECHARGE;
						bankStates[b]->stateChangeCountdown = tRP;
						break;
					case REFRESH:
					case PRECHARGE:
					case POWERDOWN_EXIT:
						bankStates[b]->currentBankState = IDLE;
						break;
					default:
						break;
				}
			}
		}
	}
	
	UpdateState();
	Step();
}

//
//Decrease countdown and send back return command to vault controller
//
void DRAM::UpdateState()
{
	if(readData != NULL) {
		dataCyclesLeft--;
		if(dataCyclesLeft == 0) {
			vaultContP->ReturnCommand(readData);
			readData = NULL;
		}
	}

	//Send back return command to vault controller
	if(readReturnCountdown.size() > 0 && readReturnCountdown[0]==0) {
		readData = readReturnDATA[0];
		dataCyclesLeft = BL;
		readReturnDATA.erase(readReturnDATA.begin());
		readReturnCountdown.erase(readReturnCountdown.begin());
	}
	
	for(int i=0; i<readReturnCountdown.size(); i++) {
		readReturnCountdown[i]--;
	}
}

//previousBankState
//Print current state in state log file
//
void DRAM::PrintState()
{
	bool printDR = false;
	for(int b=0; b<NUM_BANKS; b++) {
		if(previousBankState[b] != bankStates[b]->currentBankState) {
			printDR = true;
			break;
		}
	}
	
	if(printDR) {
		STATEN(ALI(17)<<header);
		STATEN("State");
		for(int b=0; b<NUM_BANKS; b++) {
			previousBankState[b] = bankStates[b]->currentBankState;
			STATEN("[");
			STATEN(*bankStates[b]);
			if(bankStates[b]->currentBankState == ROW_ACTIVE) {
				STATEN("-"<<bankStates[b]->openRowAddress);
			}
			STATEN("]");
		}
		STATEN(endl);
	}
}

uint64_t DRAM::get_currentclk()
{
	return currentClockCycle;
}

} //namespace CasHMC