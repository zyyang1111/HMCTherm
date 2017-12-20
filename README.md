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


## 3. Folder Directory
  - multi2sim_memtrace : architectural simulator
  - McPAT : power simulator (McPAT)
  - sources : source files for the HMC simulator 
  - result : Log files for the HMC simulation 
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
  Please see here for details of the arguments of the architectural simulation. This simulation generates four output files: 
    - **pipeline.out** : <br />
    - **mem.out** : <br />
    - **net.out** : <br />
    - **mem_trace.tr** : the trace of memory access

   - (3) generate McPAT input files using pipeline.out, mem.out and net.out
   ``` 
   $ cd $(HMCThermROOT)/CPU_power_gen/
   $ matlab test_m2s_to_mcpat.m
   ```
   This will generate mcpat.xml in $(HMCThermROOT)/Example/FFT\_example/. Please see here for more details of the usage of this MATLAB script.

   - (4) run McPAT
   ``` 
   $ cd $(HMCThermROOT)/Example/FFT_example/
   $ sh run_mcpat.sh
   ```
   This will generate mcpat.txt in $(HMCThermROOT)/Example/FFT\_example/.

   - (5) calculate the power profile for the multi-core CPU
   ```
   $ cd $(HMCThermROOT)/CPU_power_gen/
   $ matlab test_mcpattxt_to_powermap.m
   ```
   This will generate an input file to the HMC simulator with the information of the CPU layer power profile. This input file (logicP.in) will be located in $(HMCThermROOT)/Example/FFT\_example/. Please see here for more details of the usage of this MATLAB script.

   - (6) run the HMC simulator
   ``` 
   $ cd $(HMCThermROOT)/
   $ ./HMCTherm -c 400000 -t file -f ./Example/FFT_example/mem_trace.tr -x 1024 -y 1024 -q ./Example/FFT_example/logicP.in -e 20000
   ```
   In this command: <br /> 
      **-c** indicates the number of CPU cycles to be simulated <br /> 
      **-t** indicates whether the input memory trace is from a file (file) or generated randomly (random) <br />
      **-f** if the memory trace is from a file, the trace file name is specified here <br />
      **-x** and -y specifies the size (in byte) of a "mat" which is the smallest unit of memory [Default 512] <br />
      **-e** is the time step to print out the transient temperature and power profile <br /> 
      **-h** will print out the help of the HMC simulator <br /> 

   The HMC simulation will generate five output files: <br />
      - **/result/file.out** : summary of the HMC latency <br /> 
      - **temperature\_trace.csv** : the temperature traces for each time step <br /> 
      - **power\_trace.csv** : the power traces for each time step <br />
      - **Average_power.csv** : the average power profile within the time of simulation <br />
      - **static_temperature.csv** : the static thermal profile given the average power <br />
   Please see here for more details of the usage of the HMC simulator and the specifications of each file.

   - (7) Plot the power and temperature profile using the python scripts in $(HMCThermROOT)/script/
     - plot_trans.py : plot the transient power or temperature profile 
     - plot_staticT.py : plot the static temperature profile
     - plot_staticP.py : plot the average power profile



