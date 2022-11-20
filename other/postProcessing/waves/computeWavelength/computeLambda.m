%% Computes wavelength of a wave
%  Given the period, water depth and wave height
function [lambda] = computeLambda(period, depth, g)
    lambda = 0; % Variable initialization
    fun =  @(lambda)dispersionRelation(lambda, period, depth, g);
    x0 = 0; % Intial guess
    lambda = fsolve(fun,x0); % Solve the non-linear dispersion relation
end