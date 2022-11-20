%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Compute wavelength given the period and depth  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define the parameters
period = 8; % Wave period [s]
depth = 150; % Water depth [m]
g = 9.81; % Gravity acceleration [m/s^2]

%% Compute wavelength
lambda = computeLambda(period, depth, g) % Wavelength