
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Overall turbine loads as a function of time  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% Read ALM file info
case_name = "name_of_the_case\";
addpath("data\" + case_name + "turbines\0\");

filename = "turbine.csv"; % Name of the turbine

raw_data = table2array(readtable(filename));
nOuter = nnz(raw_data(:,1)==raw_data(1,1)); % Number of PIMPLE loops
time = raw_data(nOuter:nOuter:end,1); % Time vector
azimuth = raw_data(nOuter:nOuter:end,2); % Azimuth vector
% Turbine coefficients
cp = raw_data(nOuter:nOuter:end,4); % Power
cd = raw_data(nOuter:nOuter:end,5); % Drag
cq = raw_data(nOuter:nOuter:end,6); % Torque
ct = raw_data(nOuter:nOuter:end,7); % Thrust
% Blade coefficients (drag and torque)
cd1 = raw_data(nOuter:nOuter:end,11);
cq1 = raw_data(nOuter:nOuter:end,12);
cd2 = raw_data(nOuter:nOuter:end,13);
cq2 = raw_data(nOuter:nOuter:end,14);
cd3 = raw_data(nOuter:nOuter:end,15);
cq3 = raw_data(nOuter:nOuter:end,16);

rmpath("data\" + case_name + "turbines\0\");

%% Turbine variables
U = 9; %freestream
D = 126; %diameter
A = pi*(D/2)^2; % Area
dens = 1.225; %Density

drag = cd*dens*0.5*A*U^2; %Drag force in N
thrust = ct*dens*0.5*A*U^2; %Thrust force in N
torque = cq*dens*0.5*A*U^2*(D/2); %Drag force in Nm

%%  Time window to visualize
t0 = 0; % Start time
tf = 320;%time(end); % End time
[~, index0] = min(abs(time-t0));
[~, indexf] = min(abs(time-tf));

%% Mean loads
drag_mean = mean(drag(index0:indexf))
thrust_mean = mean(thrust(index0:indexf))
torque_mean = mean(torque(index0:indexf))
cd_mean = mean(cd(index0:indexf))
ct_mean = mean(ct(index0:indexf))
cp_mean = mean(cp(index0:indexf))

