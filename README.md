# HMCTemp
HMCTemp: A Cycle-accurate Simulator for Hybrid Memory Cube with thermal analysis
version HMCTemp v1.0 - 2017.06.26

1. Developer 

  Zhiyuan Yang
  Ankur Srivastava 
  University of Maryland 
  zyyang [at] umd [dot] edu

2. About HMCTemp

  HMCTemp is built based on CasHMC v1.2 (https://github.com/estwings57/CasHMC) which 
  provides cycle-accurate simulation of every module in HMC and generates analysis 
  results. HMCTemp adds a power-thermal module to CasHMC such that the power profile 
  can be captured. After the simulation, HMCTemp can generate the steady-state 
  temperature profile of the HMC. Given the time step, HMCTemp can provide transient 
  thermal simulation. 

  HMCTemp is basically implemented in C++. The integrated grid-based thermal solver is 
  based on SuperLU that significantly speeds up the thermal simulation. HMCTemp can take
  the power profile of the processor which is stacked with the HMC into the thermal 
  simulation. The processor power profile is specified in logicP.in

3. Folder directory

  graph : Gnuplot script file and a graph data file after a simulation run is over
  result : Log files after a simulation run is over
  sources : All CasHMC source files
  trace : The example of trace files (The files are extracted from SPEC CPU2006 benchmarks)
  power_trace: power profile for each time step
  temperature_trace: temperature profile for each time step 

4. Build HMCTemp 

  (1) Dependency. 
  The integrated grid-based thermal solver is based on SuperLU_MT (SuperLU for shared memory
  parallel machines). SuperLU_MT should be installed before installing HMCTemp. SuperLU_MT 
  can be downloaded from http://crd-legacy.lbl.gov/~xiaoye/SuperLU/

  (2) Modify Makefile
  In Makefile, modify the SuperLUroot to the directory where SuperLU_MT is installed in your
  machine

  (3) Install 
  To build an HMCTemp
  $ make

5. Running HMCTemp
  
  (1) Command line arguments
  -c (--cycle)   : The number of CPU cycles to be simulated
  -t (--trace)   : Trace type ('random' or 'file')
  -u (--util)    : Requests frequency (0 = no requests, 1 = as fast as possible) [Default 0.1]
  -r (--rwratio) : (%) The percentage of reads in request stream [Default 80]
  -f (--file)    : Trace file name
  -h (--help)    : Simulation option help
  -x (--gridx)   : number of grids in a bank in x direction [Default 1]
  -y (--gridy)   : number of grids in a bank in y direction [Default 1]

  (2) Example
  $ ./HMCTemp -c 100000 -t random -u 0.1 -r 60 -x 4 -y 4
