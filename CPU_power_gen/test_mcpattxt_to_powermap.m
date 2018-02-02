close all 
clear
clc

addpath ./TCG_operations/;
addpath ./matlab_bgl/;
addpath ./lib;

dir_name = ['../Example/FFT_example/'];
siminfo_dir_name = ['../Example/sim_info/'];

% %%%%%%%%%%%%%%%%%%%%% Predefined Parameters %%%%%%%%%%%%%%%%%%%%
name_vector={ ...
    'bpred', ...
    'BTB', ...
    'I$', ...
    'fetchQ', ...
    'decode', ...
    'rename', ...
    'ROB', ...
    'IQ', ...
    'LSQ', ...
    'RF', ...
    'D$', ...
    'EX', ...
    'router', ...
    'L2', ...
    'MMU', ...
    'MC' ...
    };
section_vector={ ...
    {'Branch Predictor'}, ...
    {'Branch Target Buffer'}, ...
    {'Instruction Cache'}, ...
    {'Instruction Buffer'}, ...
    {'Instruction Decoder'}, ...
    {'Renaming Unit'}, ...
    {'ROB'}, ...
    {'Instruction Window','FP Instruction Window'}, ...
    {'LoadQ', 'StoreQ'}, ...
    {'Register Files'}, ...
    {'Data Cache'}, ...
    {'Integer ALUs (Count','Floating Point Units (FPUs) (Count','Complex ALUs (Mul/Div) (Count'}, ...
    {'Router'}, ...
    {'L2'}, ...
    {'Memory Management Unit'}, ...
    {'Memory Controller'} ...
    };

% num_cores --> average number
% cores
a_cores = 16;
num_entries_vector=[...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    a_cores, ...
    1, ...
    a_cores, ...
    a_cores, ...
    1 ...
    ];

core_per_MC = 4; % number of cores per memory controller 
num_cores = 16; % number of cores 
num_MC = num_cores / core_per_MC;
num_mod=length(name_vector);
aspect=rand(num_mod,1)*4+1;

freq = 2; % GHz -- operating frequency 
tech=32; % nm -- fab technology node
vdd_real=0.9; % V -- supply voltage 
wire_level=2; %1: M1, 2: intermediate 3: global
wire_delay=calc_wire_delay( tech, vdd_real, wire_level ); %ps/mm
single_cycle_WL=1e3/freq/wire_delay; %milimeters

% %%%%%%%%%%%%%%%%%%%%%%%%% Main Calculation %%%%%%%%%%%%%%%%%%%%%%%

[ section section_name ] = read_mcpat_power_area( [dir_name 'mcpat.txt'] );

for i=1:length(section_vector)
    area(i,1)=0;
    dyn(i,1)=0;
    leak(i,1)=0;
    for j=1:length(section_vector{i})
        area(i)=area(i)+str2num(reference_multi(section,section_name,section_vector{i}{j},'Area',1))*1e6; %um^2
        for k=1:num_entries_vector(i)
            dyn(i)=dyn(i)+str2num(reference_multi(section,section_name,section_vector{i}{j},'Runtime Dynamic',k));
            leak(i)=leak(i)+str2num(reference_multi(section,section_name,section_vector{i}{j},'Subthreshold Leakage',k));
            leak(i)=leak(i)+str2num(reference_multi(section,section_name,section_vector{i}{j},'Gate Leakage',k));
        end
    end
    dyn(i)=dyn(i)/num_entries_vector(i);
    leak(i)=leak(i)/num_entries_vector(i);
end

area(find(strcmp(name_vector,'MC')))=area(find(strcmp(name_vector,'MC')))/core_per_MC;
dyn(find(strcmp(name_vector,'MC')))=dyn(find(strcmp(name_vector,'MC')))/core_per_MC;
leak(find(strcmp(name_vector,'MC')))=leak(find(strcmp(name_vector,'MC')))/core_per_MC;
                                            
power_den=(dyn+leak)./area.*1e8; %W/cm^2

width=sqrt(area.*aspect);
hight=width./aspect;

load([siminfo_dir_name 'manual_FP.mat']);
[ x, y ] = calc_FP( Gh, Gv, width, hight );

% ========= export the power information for each core ============
csv_data=[width hight x y dyn leak ones(size(x))];
dlmwrite([dir_name 'initial_floorplan.csv'],csv_data,'precision',10);

% ====== generate the power density map for the whole circuit ======
% use the precalculated floorplan infomation 
csv_data_dim=csvread([siminfo_dir_name 'final_floorplan_reliability_unaware_aircooled_modsim.csv']);
width=csv_data_dim(:,1);
height=csv_data_dim(:,2);
x=csv_data_dim(:,3);
y=csv_data_dim(:,4);
P_leak(:)=csv_data_dim(:,6);

P=csv_data(:,5);
L=csv_data(:,7:end)';

P=mean(P,2);
P_leak=mean(P_leak,2);
P_den=P./(width.*height*1e-12);
P_leak_den=P_leak./(width.*height*1e-12);
area=width.*height;
aspect=width./height; % because we use the precalculated parameter, aspect ratio should be recalculated

MC_idx=find(strcmp(name_vector,'MC'));
core_idx=setdiff(1:length(name_vector),MC_idx);
MC_area=sum(width(MC_idx).*height(MC_idx))*1e-6; %convert um^2 to mm^2
core_area=sum(width(core_idx).*height(core_idx))*1e-6; %convert um^2 to mm^2
                                            
topologies = enumerate_network_topologies( num_cores );
[ topology, topology_delay ] = identify_best_topology( topologies, num_MC, core_area, MC_area, single_cycle_WL );

% Gh and Gv are from "manual_FP.mat"
[ x, y ] = calc_FP( Gh, Gv, width, height );
dimX0=max(x+width)-min(x);
dimY0=max(y+height)-min(y);
                                            
diag0=sqrt((dimX0*topology(2)/topology(1))^2+dimY0^2);
A0=dimX0*dimY0;

A_opt=sum(width.*height);

diag_opt=sqrt(A_opt*topology(2)/topology(1));

A_norm0=A0/A_opt;
diag_norm0=diag0/diag_opt;

x=round(x/10)*10;
y=round(y/10)*10;
width=round(width/10)*10;
height=round(height/10)*10;

[ T0, powerden_dyn, powerden_leak ] = generate_Pden_map( x, y, width, height, P_den, P_leak_den, topology);
ChipX = max(x + width) * 4;
ChipY = max(y + height) * 4; 

powerden = powerden_dyn + powerden_leak;
gridx = ChipX / size(powerden, 1) * 1e-6; % [m]
gridy = ChipY / size(powerden, 2) * 1e-6; % [m]

x_size = size(powerden,1) / 15; 
y_size = size(powerden,1) / 15; 
gridx_n = ChipX / x_size * 1e-6; 
gridy_n = ChipY / y_size * 1e-6; 
powerden_n = imresize(powerden, [x_size, y_size]); 

powerM = powerden_n * gridx_n * gridy_n; 

% imagesc(powerM)

fw = fopen([dir_name, 'logicP.in'], 'w'); 
fprintf(fw, '%d %d\n', size(powerM,1), size(powerM,2)); 
for i = 1 : size(powerM,1)
    for j = 1 : size(powerM, 2)
        fprintf(fw, '%.7f ', powerM(i,j)); 
    end
end
fprintf(fw, '\n');
fclose(fw);

% sum(powerM(:))