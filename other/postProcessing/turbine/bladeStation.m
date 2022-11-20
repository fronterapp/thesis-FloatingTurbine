
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Get variables at one blade station as a function of time  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

case_name = "name_of_the_case\";
addpath("data\" + case_name + "actuatorLineElements\0\");

%% Read ALM file info
blades = [1 2 3]; % Search for blade1, blade2 and blade3
element_num = "0"; %Actuator element defining the radial position of interest

for i=1:length(blades)
    filename = "turbine.blade" + int2str(blades(i)) + ".element" + element_num + ".csv";
    raw_data(:,:,i) = table2array(readtable(filename));
end
rmpath("data\" + case_name + "actuatorLineElements\0\"); % Remove path

%% Time vector
time_all = raw_data(:,1,1); % Time vector
nOuter = nnz(time_all==time_all(1)); % Number of PIMPLE loops
time = time_all(nOuter:nOuter:end); %Time vector, last PIMPLE loop only

%% Access variables
root_d = raw_data(1,2,1); % Vector of radial positions

for i=1:length(blades)
    % Position
    x(:,i) = raw_data(nOuter:nOuter:end,3,i);
    y(:,i) = raw_data(nOuter:nOuter:end,4,i);
    z(:,i) = raw_data(nOuter:nOuter:end,5,i);
    % Velocity magnitude
    vel(:,i) = raw_data(nOuter:nOuter:end,6,i);
    % Prescribed motion velocity
    pres_vx(:,i) = raw_data(nOuter:nOuter:end,7,i);
    pres_vy(:,i) = raw_data(nOuter:nOuter:end,8,i);
    pres_vz(:,i) = raw_data(nOuter:nOuter:end,9,i);
    % Relative velocity
    rel_vx(:,i) = raw_data(nOuter:nOuter:end,23,i);
    rel_vy(:,i) = raw_data(nOuter:nOuter:end,24,i);
    rel_vz(:,i) = raw_data(nOuter:nOuter:end,25,i);
    rel_mag(:,i) = sqrt(rel_vx(:,i).^2 + rel_vy(:,i).^2 + rel_vz(:,i).^2);
    % Inflow velocity
    ux(:,i) = raw_data(nOuter:nOuter:end,26,i);
    uy(:,i) = raw_data(nOuter:nOuter:end,27,i);
    uz(:,i) = raw_data(nOuter:nOuter:end,28,i);
    inflow(:,i) = sqrt(ux(:,i).^2 + uy(:,i).^2 + uz(:,i).^2);
    % Blade rotation velocity
    rot_x(:,i) = raw_data(nOuter:nOuter:end,29,i);
    rot_y(:,i) = raw_data(nOuter:nOuter:end,30,i);
    rot_z(:,i) = raw_data(nOuter:nOuter:end,31,i);
    rot_mag(:,i) = sqrt(rot_x(:,i).^2 + rot_y(:,i).^2 + rot_z(:,i).^2);
    %Reynolds, angle of attack, aero. coeffs
    re(:,i) = raw_data(nOuter:nOuter:end,10,i);
    alpha(:,i) = raw_data(nOuter:nOuter:end,11,i);
    alpha_geo(:,i) = raw_data(nOuter:nOuter:end,12,i);
    cl(:,i) = raw_data(nOuter:nOuter:end,13,i);
    cd(:,i) = raw_data(nOuter:nOuter:end,14,i);
end
