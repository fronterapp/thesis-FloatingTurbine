%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  Get time-averaged variables along the span  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

case_name = "name_of_the_case\";
addpath("data\" + case_name + "actuatorLineElements\0\");

%% Read ALM file info
U_inf = 4.19; % Freestream vel
blades = [1 2 3]; % Blades we want to sample
n_el = 32; %Number of actuator elements

%% Acess time vector
filename_dummy = "turbine.blade1.element1.csv";%File from which to get the time vector
raw_data_0(:,:) = table2array(readtable(filename_dummy));
time = raw_data_0(:,1);

%% Compute average
t0 = 0; % Start time
tf = 1; % End time
[~, index0] = min(abs(time-t0));
[~, indexf] = min(abs(time-tf));

%% Access files and average
for i=1:length(blades)
    for j=1:n_el
        filename = "turbine.blade" + int2str(blades(i)) + ".element" + int2str(j-1) + ".csv";
        data_array = table2array(readtable(filename));
        all_raw_data(j,:,:,i) = data_array(index0:indexf,:);
        all_raw_data2(j,:,:,i) = data_array;
        %Now average in time
        raw_data(j,:,i) = mean(all_raw_data(j,:,:,i));       
    end
end

%% Compute rotor-average induction factor
ref_r = raw_data(:,2,1); % Radial positions
ref_areas = pi*(ref_r(2:end).^2-ref_r(1:end-1).^2); % Areas
for i=1:length(blades)
    for j=1:length(time)
        ux_radial = all_raw_data2(:,j,26,i);
        a_radial = (U_inf - ux_radial ) / U_inf;
        mean_a_radial = 0.5*(a_radial(2:end)+a_radial(1:end-1));
        induction(j,i) = sum(mean_a_radial .* ref_areas)./sum(ref_areas);
    end
end

rotor_a = mean(induction,2);

%% Access averaged variables
for i=1:length(blades)
    % Position
    r(:,i) = raw_data(:,2,i);
    x(:,i) = raw_data(:,3,i);
    y(:,i) = raw_data(:,4,i);
    z(:,i) = raw_data(:,5,i);
    % Velocity magnitude
    vel(:,i) = raw_data(:,6,i);
    % Prescribed motion velocity
    pres_vx(:,i) = raw_data(:,7,i);
    pres_vy(:,i) = raw_data(:,8,i);
    pres_vz(:,i) = raw_data(:,9,i);
    % Relative velocity
    rel_vx(:,i) = raw_data(:,23,i);
    rel_vy(:,i) = raw_data(:,24,i);
    rel_vz(:,i) = raw_data(:,25,i);
    rel_mag(:,i) = sqrt(rel_vx(:,i).^2 + rel_vy(:,i).^2 + rel_vz(:,i).^2);
    % Inflow velocity
    ux(:,i) = raw_data(:,26,i);
    uy(:,i) = raw_data(:,27,i);
    uz(:,i) = raw_data(:,28,i);
    inflow(:,i) = sqrt(ux(:,i).^2 + uy(:,i).^2 + uz(:,i).^2);
    a(:,i) = (U_inf - ux(:,i) ) / U_inf; %Axial induction
    % Blade rotation velocity
    rot_x(:,i) = raw_data(:,29,i);
    rot_y(:,i) = raw_data(:,30,i);
    rot_z(:,i) = raw_data(:,31,i);
    rot_mag(:,i) = sqrt(rot_x(:,i).^2 + rot_y(:,i).^2 + rot_z(:,i).^2);
    %Reynolds, angle of attack, aero. coeffs
    re(:,i) = raw_data(:,10,i);
    alpha(:,i) = raw_data(:,11,i);
    alpha_geo(:,i) = raw_data(:,12,i);
    cl(:,i) = raw_data(:,13,i);
    cd(:,i) = raw_data(:,14,i);
end