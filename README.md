# HMCTherm
HMCTherm: A Cycle-accurate Simulator for Hybrid Memory Cube with thermal analysis
version HMCTherm v1.0 - 2017.06.26

## 1. Developer 

  Zhiyuan Yang <br />
  Ankur Srivastava <br />
  University of Maryland <br />
  zyyang [at] umd [dot] edu

## 2. About HMCTherm

  HMCTherm is a comprehensive simulation framework for Stacked-Memory-on-CPU architecture. Given the architectural description of multi-core CPUs and HMCs, HMCTherm can simulate the 3D thermal profile (both transient and static) of the HMC when a certain program is running. 

  HMCTherm is composed of: 
  - An HMC (Hybrid-Memory-Cube) Simulator based on CasHMC v1.2 (https://github.com/estwings57/CasHMC) which takes memory traces as inputs and performs cycle-accurate simulation of HMC to get the delay, power and thermal information. 
  - An architectural simulator for multi-core CPUs based on Multi2Sim (https://github.com/Multi2Sim/multi2sim) which performs architectural simulation and generates memory trace.
  - A power simulator based on McPAT (https://github.com/HewlettPackard/mcpat)
  - Other Shell script files which fulfill the interface between different simulators 

  A grid-based thermal solver is integrated to HMCTherm which is based on SuperLU\_MT (http://crd-legacy.lbl.gov/~xiaoye/SuperLU/). SuperLU is a general purpose library for the direct solution of large, sparse, nonsysmmetric systems of linear equations; and SuperLU\_MT is a version that enables multi-threaded computation. 


3. Folder directory

  graph : Gnuplot script file and a graph data file after a simulation run is over
  result : Log files after a simulation run is over
  sources : All HMCTherm source files
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
  To build an HMCTherm
  $ make

5. Running HMCTherm
  
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
  $ ./HMCTherm -c 100000 -t random -u 0.1 -r 60 -x 4 -y 4
