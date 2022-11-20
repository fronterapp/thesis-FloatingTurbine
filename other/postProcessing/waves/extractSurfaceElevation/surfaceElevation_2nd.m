%% Compute surface elevation according to Stoke's second-order theory
%  Inputs are:
%  Position (x), time (t), wavelength(lambda)
%  Period (T), wave height (H), and water depth (h)
%  Output is surface elevation (eta)
function [eta] = surfaceElevation_2nd(x, t, lambda,T, H, h)
    a = H/2; % Wave amplitude
    k  = 2*pi/lambda; % Wavenumber
    w = 2*pi/T; %Wave frequency
    theta = k*x - w*t; % Phase
    sigma = tanh(k*h);
    f = (3-sigma^2)/(4*sigma^3);
    eta = a*(cos(theta)+ k*a*cos(2*theta)*f); % Surface elevation
end