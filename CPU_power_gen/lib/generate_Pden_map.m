function [powerden_dyn, powerden_leak, powerden_arbitrary_dyn ] = generate_Pden_map( x, y, w, h, P, P_leak, topology)

%     dram_power=5;
    
    %%%%x,y,w,h are in units of microns, convert to units of 10-microns
    x=round(x/10);
    y=round(y/10);
    w=round(w/10);
    h=round(h/10);

    powerden_dyn=[];
    powerden_leak=[];
    
%     n=length(x);
% %     k=size(L,1);
% 
%     num_DRAM_layers=4; %number of transistor layers for DRAM
% %     num_logic_layers=num_cores/16; %number of transistor layers for logic
%     num_logic_layers=topology(3);
%     
% %     v=0.5; %fluid velocity
%     FR=1e-7*num_logic_layers; %this means v=0.5 in initial 32/64 core FP with all channels on
% %     pumping_power=10e-3*num_logic_layers; 
%     pumping_power=power_per_layer*num_logic_layers; 
%     thresh=1; %stopping criteria: maximum difference in T between iterations
%     max_iter=10; %maximum iterations of thermal-leakage loop
% %     cooled=1;
% 
%     Tamb = 40; %ambient temp in deg C
%     Tin = 25; %fludic input temp in deg C
% 
%     transient=0;
%     tspan=[];
%     thermal_runaway=0;
% 
%     %%%%LEAKAGE MODEL
%     leakage_model_c2=5.121;
%     leakage_model_c1=-6.013;
%     leakage_model_c0=1.892;

%     %**************************************MATERIAL PROPERTIES****************
% 
%     %%%%WATER%%%
%     Kf = 0.6069;                    % thermal conductivity of coolant
%     ru = 997;                       % density of coolant, unit: kg/m^3
%     Cp = 4181.7;                    % specific heat of coolant
%     Cmc = 4.1796e6;                 % volumetric heat capacity
%     viscosity = 8.899e-4;           % viscosity of coolant
% 
%     %%%%SILICON%%%
%     Ksi=148;
%     Csi=1.66*1e6;
% 
%     %%%%OXIDE%%%
%     Kox=1.4;
%     Cox=1.65*1e6;
% 
%     %%%%WIREING%%%
%     Kwire=2.25; %this number is from 3D ICE for BEOL
%     Cwire=2*1e6; %I made this up
% 
%     %%%%COPPER%%%
%     Kcu=401;
%     Ccu=3.2*1e6;
% 
%     %**************************************STACK PARAMETERS****************
% 
%     %layer heights in meters
%     Hmc_dummy=20e-6;
%     if cooled==1
%         Hmc=200e-6;
%     else
%         Hmc=Hmc_dummy;
%     end
%     Hwire=15e-6;
%     Hact=5e-6;
%     Hbase=55e-6;
%     Hox=15e-6;
%     Hwafer=1000e-6-Hact-Hbase-Hmc_dummy;
% 
%     %**************************************FLOORPLAN DESCRIPTION****************

%     Wmc=100e-6;
%     Wmc=200e-6;

%     gridX=Wmc;
%     gridZ=1000e-6;
%     gridZ=800e-6;
%     gridZ=600e-6;
%     gridZ=400e-6;
%     gridZ=200e-6;
    
    % the values are set by yzy 5/18/2016
    gridX = 10e-6;
    gridZ = 10e-6;

%     rTSV = 2.5e-6;
%     Ktsv = Kcu;
% 
%     %**************************************CONSTANTS****************
%     Dh = 4 * Wmc * Hmc / 2 / (Wmc + Hmc);   % hydrolic diameter
%     ar = Hmc / Wmc;                         % aspect ratio
%     const = (4.7 + 19.64 * ((ar) ^ 2 + 1) / (ar + 1) ^ 2); %this constant is refered to as "gamma" in publications
%     correct = 5; %correction constant, do not change
%     T0 = 273; %0 deg C in kelvins
% 
%     %**************************************CREATE STACK****************

%     H=[];
%     K=[];
%     C=[];
%     layerP=[];
%     tierMap=[];
%     mapTSV=[];
% 
%     num_DRAM_layers=4; %number of transistor layers for DRAM
% %     num_logic_layers=num_cores/16; %number of transistor layers for logic
    num_logic_layers=topology(3);

    %layer vectors start from top (heatsink) to bottom (pins)
%     for layer_num=1:num_DRAM_layers
%         if layer_num==1                 %top layer
%             tierMap=[tierMap layer_num*ones(1,6)];
%             layerP=[layerP length(H)+4];
%             mapTSV=[mapTSV 0*ones(1,6)];
%             H=[H Hwafer Hmc_dummy Hbase Hact Hwire Hox];
%             C=[C Csi Csi Csi Csi Cwire Cox];
%             K=[K Ksi Ksi Ksi Ksi Kwire Kox];
%         else                            %intermediate layer
%             tierMap=[tierMap layer_num*ones(1,5)];
%             layerP=[layerP length(H)+3];
%             mapTSV=[mapTSV 1*ones(1,5)];
%             H=[H Hmc_dummy Hbase Hact Hwire Hox];
%             C=[C Csi Csi Csi Cwire Cox];
%             K=[K Ksi Ksi Ksi Kwire Kox];
%         end
%     end
% 
%     for layer_num=1:num_logic_layers
%         if layer_num==num_logic_layers  %bottom layer
%             tierMap=[tierMap layer_num+num_DRAM_layers*ones(1,4)];
%             layerP=[layerP length(H)+3];
%             mapTSV=[mapTSV 1*ones(1,4)];
%             H=[H Hmc Hbase Hact Hwire];
%             C=[C Csi Csi Csi Cwire];
%             K=[K Ksi Ksi Ksi Kwire];
%         else                            %intermediate layer
%             tierMap=[tierMap layer_num+num_DRAM_layers*ones(1,5)];
%             layerP=[layerP length(H)+3];
%             mapTSV=[mapTSV 1*ones(1,5)];
%             H=[H Hmc Hbase Hact Hwire Hox];
%             C=[C Csi Csi Csi Cwire Cox];
%             K=[K Ksi Ksi Ksi Kwire Kox];
%         end
%     end
% 
%     numLayer=length(H); %number of material layers
%     numP=num_DRAM_layers+num_logic_layers;
    numP = num_logic_layers; 
%     numMC=numP; %this is a requirement due to coding error
%     layerMC=layerP-2; %each MC layer is 2 above active layer (base layer in between)

    %**************************************CREATE POWER MAP****************
    tile_num_rows=max(x+w)-min(x);
    tile_num_cols=max(y+h)-min(y);

%     powerden_arbitrary_dyn=zeros(tile_num_rows,tile_num_cols,num_DRAM_layers+num_logic_layers);
%     powerden_arbitrary_leak=zeros(tile_num_rows,tile_num_cols,num_DRAM_layers+num_logic_layers);
    
    powerden_arbitrary_dyn=zeros(tile_num_rows,tile_num_cols,num_logic_layers);
    powerden_arbitrary_leak=zeros(tile_num_rows,tile_num_cols,num_logic_layers);
    
    for i=1:num_logic_layers
        for j=1:length(x)
            powerden_arbitrary_dyn(x(j)+(1:w(j)),y(j)+(1:h(j)),i)=P(j);
            powerden_arbitrary_leak(x(j)+(1:w(j)),y(j)+(1:h(j)),i)=P_leak(j);
        end
    end
    
    %%%%tile the powerden
    powerden_arbitrary_dyn=repmat(powerden_arbitrary_dyn,[topology(1),topology(2),1]);
    powerden_arbitrary_leak=repmat(powerden_arbitrary_leak,[topology(1),topology(2),1]);
    
    x_pdn_step = size(powerden_arbitrary_dyn, 1) / 4;
    y_pdn_step = size(powerden_arbitrary_dyn, 2) / 4;
    
    for i = 1 : 4
        for j = 1 : 4
            x_span = (i-1)*x_pdn_step + 1 : i*x_pdn_step;
            y_span = (j-1)*y_pdn_step + 1 : j*y_pdn_step; 
            ind = (i-1)*4 + j;
        end
    end

    
    %calculate true dynamic and leakage power
    true_dyn=sum(sum(sum(powerden_arbitrary_dyn*10e-6*10e-6)));
    true_leak=sum(sum(sum(powerden_arbitrary_leak*10e-6*10e-6)));
    
    
    %calculate chip size
    chip_width=size(powerden_arbitrary_dyn,1)/1e5; %convert to meters
    chip_height=size(powerden_arbitrary_dyn,2)/1e5;
    
%     %%%%%%%FLIP CHIP VERTICALLY
%     figure
%     imagesc(powerden_arbitrary_dyn(:,:,1))
%     powerden_arbitrary_dyn=flip(powerden_arbitrary_dyn,3);
%     powerden_arbitrary_leak=flip(powerden_arbitrary_leak,3);
%     figure
%     imagesc(powerden_arbitrary_dyn(:,:,1))

    %**************************************CREATE FLOORPLAN****************
    %normalize chip_width and chip_height to nearest grid size
    chip_width=round(chip_width/gridX)*gridX;
    chip_height=round(chip_height/gridZ)*gridZ;


    %this is by definition
    dimX=fix(chip_width/gridX);
    dimZ=fix(chip_height/gridZ);

    %********************************resize powerden_arbitrary to powerden (dimX by dimZ)******************
    clear powerden_dyn powerden_leak
    for layer_num=1:numP
        powerden_dyn(:,:,layer_num)=imresize(powerden_arbitrary_dyn(:,:,layer_num),[dimX dimZ],'bilinear');
        powerden_leak(:,:,layer_num)=imresize(powerden_arbitrary_leak(:,:,layer_num),[dimX dimZ],'bilinear');
    end
    
    %calculate dynamic and leakage power after resizeing
    resized_dyn=sum(sum(sum(powerden_dyn*gridX*gridZ)));
    resized_leak=sum(sum(sum(powerden_leak*gridX*gridZ)));

    %scale dynamic and/or leakage power to correct resizing error
    powerden_dyn=powerden_dyn*true_dyn/resized_dyn;
%     powerden_leak=powerden_dyn*true_leak/resized_leak;



end

