%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read free-surface elevation data and compute error with Stokes theory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% READ DATA
filename = "SE_2nd_order.dat";
addpath('wave_data')

raw_data = table2array(readtable(filename));

idx=find(raw_data(:,1)>0); %Gives indices of positive values of time
data = raw_data(idx,:);
time = data(:,1); % Time vector
surf_elevation = data(:,2:end); % Surface elevation

%% WAVE PARAMETERS
wave_period = 8; % Wave period [s]
height_ratio = 0.005; % Steepness, wave height to wavelength ratio
depth_ratio = 0.4; %Water depth to length ratio
wave_length = 9.81/(2*pi)*wave_period^2*tanh(2*pi*depth_ratio); % Dispersion relation

wave_height = height_ratio*wave_length; % Wave height [m]
water_depth = depth_ratio*wave_length; % Water depth [m]
wave_freq = 2*pi/wave_period; %Wave frequency [rad/s]

%% CONSTRUCT REFERENCE WAVE
n_periods = time(end)/wave_period; % Number of periods in the data
n_samples = length(time)/n_periods; % Mean number of samples per period
ref_time = linspace(0,wave_period,n_samples); % Reference time for cosine wave
gauges_pos = linspace(0,2*wave_length,size(surf_elevation,2)); %Dimensionless position of the gauges (x/lambda)

% Compute second-order phase average
ref_phase_avg = surfaceElevation_2nd(0, ref_time, wave_length, wave_period, wave_height, water_depth);

% Compute the Stokes second-order surface elevation at evry position and time
for i=1:length(time)
    ref_wave_matrix(i,:) = surfaceElevation_2nd(gauges_pos, time(i), wave_length, wave_period, wave_height, water_depth);
end

%% COMPUTE ERROR ALONG THE DOMAIN
error_amp = zeros(size(surf_elevation,2),1);
error_harm = zeros(size(surf_elevation,2),1);

for i=1:size(surf_elevation,2)
     [phase_avg(:,i)] = averagePhase(surf_elevation(:,i),0,n_samples);
     % nRMSE error, normaslised with wave height
     nrmse(i) = (100/wave_height)*sqrt(mean((phase_avg(:,i)' - ref_phase_avg).^2)); 
     % relative error, normaslised with wave height
     rel_error(:,i) = (100/wave_height)*(surf_elevation(:,i) - ref_wave_matrix(:,i));
     % relative error of the phase average, normaslised with wave height
     avg_rel_error(:,i) = (100/wave_height)*(phase_avg(:,i)' - ref_phase_avg);
end
