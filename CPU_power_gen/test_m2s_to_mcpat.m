close all
clear
clc

addpath ./lib;

% This script is for testing "test_m2s_to_mcpat"
dir_name = ['../Example/FFT_example/'];
siminfo_dir_name = ['../Example/sim_info/'];

mem_file=[dir_name 'mem.out'];
pipeline_file=[dir_name 'pipeline.out'];
memconfig_file=[siminfo_dir_name 'memconfig'];
net_file=[dir_name 'net-l1-l2_net.out'];
xml_file=[dir_name 'mcpat.xml'];

NOC_width= 128; % Byte -- NOC bandwidth
core_per_MC = 4; % number of cores per memory controller
num_cores = 16; % number of cores -- should be the same as the configuration file 
num_MC=num_cores/core_per_MC;
clk=2000; % MHz -- CPU frequency
mem_clk = 1250; % MHz -- memory frequency 
mem_size = 4096; % MB -- HMC size 
num_banks = 16; % number of banks 
network_latency = 3; 
router_ports = 5;

m2s_to_mcpat( mem_file, pipeline_file, memconfig_file, net_file, xml_file, num_MC, clk, router_ports, network_latency, NOC_width, mem_clk, mem_size, num_banks );