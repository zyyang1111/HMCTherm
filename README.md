# HMCTherm
HMCTherm: A Cycle-accurate Simulator for Hybrid Memory Cube with thermal analysis
version HMCTherm v1.0 - 2018.04.25

## 1. Developer 

  Zhiyuan Yang <br />
  Ankur Srivastava <br />
  University of Maryland <br />
  zyyang [at] umd [dot] edu

## 2. About HMCTherm

  HMCTherm is a comprehensive simulation framework for Stacked-Memory-on-CPU architecture. Given the architectural description of multi-core CPUs and HMCs, HMCTherm can simulate the 3D thermal profile (both transient and static) of the HMC when a certain program is running. The HMCTherm is built based on CasHMC v1.2 (https://github.com/estwings57/CasHMC) which takes memory traces as inputs and performs cycle-accurate simulation of HMC to get the delay, power and thermal information. 

  A grid-based thermal solver is integrated to HMCTherm which is based on SuperLU\_MT (http://crd-legacy.lbl.gov/~xiaoye/SuperLU/). SuperLU is a general purpose library for the direct solution of large, sparse, nonsysmmetric systems of linear equations; and SuperLU\_MT is a version that enables multi-threaded computation. 


## 3. Folder Directory
  - sources : source files for the HMC simulator 
  - CPU\_power\_gen : MATLAB scripts generating interface scripts between the architectural simulator, HMC simulator and McPAT
  - Example : containing benchmarks (splash2) and scripts for running the benchmarks
  - script : containing python scripts for ploting final power and thermal profiles

## 4. Usage 

### Dependency
  The integrated grid-based thermal solver is based on SuperLU\_MT. SuperLU\_MT should be installed before installing HMCTemp. SuperLU\_MT can be downloaded from http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ and following the instruction to install the library. 

### Installation
  - (1) install the architectural simulator 
  ```
  $ cd multi2sim_memtrace
  $ make
  ```
  - (2) install McPAT
  ```
  $ cd McPAT
  $ make
  ```
  - (3) in the Makefile for the HMC simulator (under the root directory of HMCTherm), modify the SuperLUroot to where the SuperLU is installed in your local system 

  - (4) install the HMC simulator 
  ``` 
  $ make 
  ```
### Run your Simulation
  In this branch, we integrate a run-time refresh control module to the HMC and develop a simulation framework to fulfill t he simulation. This simulation framework involves Multi2Sim for function simulation of multi-core processors, McPAT for estimation of the the processor power, HMCTherm for HMC simulation and MATLAB scripts for processing interface data. 

  In order to illustrate how HMCTherm works, we run FFT (./Example/splash2/codes/kernels/fft/) and plot the power and thermal profile of HMC. 
  - (1) build FFT
  ```
  $ cd $(HMCThermROOT)/Example/splash2/codes/kernels/fft/
  $ make
  ```
  - (2) run the architectural simulation
  ```  
  $ cd $(HMCThermROOT)/Example/FFT_example/
  $ sh run.sh
  ```
  We modify the original Multi2Sim such that it can output memory trace. You can find and download this program from https://github.com/zyyang1111/multi2sim-memtrace. Be sure to change use the correct Multi2Sim directory in run.sh. The simulation requires three input files: <br />
    - **cpuconfig** <br />
    - **memconfig** <br />
    - **netconfig** <br />
We have included exmaple input files in $(HMCThermROOT)/Example/sim_info/. The users can create their own files. This simulation generates four output files: <br />
    - **pipeline.out** <br />
    - **mem.out** <br />
    - **net.out** <br />
    - **mem_trace.tr** : the trace of memory access

   - (3) generate McPAT input files using pipeline.out, mem.out and net.out
   ``` 
   $ cd $(HMCThermROOT)/CPU_power_gen/
   $ matlab test_m2s_to_mcpat.m
   ```
   This will generate mcpat.xml in $(HMCThermROOT)/Example/FFT\_example/. Please read the comments in the script for details of the variable definition. The users can use their own method to generate the input file (i.e. mcpat.xml) of McPAT.

   - (4) run McPAT
   ``` 
   $ cd $(HMCThermROOT)/Example/FFT_example/
   $ sh run_mcpat.sh
   ```
   This will generate mcpat.txt in $(HMCThermROOT)/Example/FFT\_example/. You can download McPAT from https://github.com/HewlettPackard/mcpat. When using run_mcpat.sh, be sure to change the directory of McPAT. 

   - (5) calculate the power profile for the multi-core CPU
   ```
   $ cd $(HMCThermROOT)/CPU_power_gen/
   $ matlab test_mcpattxt_to_powermap.m
   ```
   This will generate an input file to the HMC simulator with the information of the CPU layer power profile. This input file (logicP.in) will be located in $(HMCThermROOT)/Example/FFT\_example/. This power file of the CPU layer will be used for HMC simulation. lease read the comments in the script for details of the variable definition. The users can use their own method to generate the CPU power profile. 

   - (6) run the HMC simulator
   ``` 
   $ cd $(HMCThermROOT)/
   $ ./HMCTherm -c 4000000 -t file -f ./Example/FFT_example/mem_trace.tr -x 1024 -y 1024 -e 200000 -d ./Example/FFT_example -q ./Example/FFT_example/logicP.in -s ./DRAM_RT.txt -j 1
   ```
   In this command: <br /> 
      **-c** indicates the number of CPU cycles to be simulated <br /> 
      **-t** indicates whether the input memory trace is from a file (file) or generated randomly (random) <br />
      **-f** if the memory trace is from a file, the trace file name is specified here <br />
      **-x** and **-y** specifies the size (in byte) of a "mat" which is the smallest unit of memory [Default 512] <br />
      **-e** is the time step to print out the transient temperature and power profile <br /> 
      **-d** specifies the directory saving the output files
      **-q** specify the file of power profile of the processor layer
      **-j** controls whether the run-time control is performed [1] or not [0]. [Default 0]
      **-h** will print out the help of the HMC simulator <br /> 


