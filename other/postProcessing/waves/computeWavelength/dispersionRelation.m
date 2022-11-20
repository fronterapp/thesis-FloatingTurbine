%% Wave dispersion relation for shallow waters
%  Returns zero if the dispersion relation is met
function [F] = dispersionRelation(lambda, period, depth, g)
    F = lambda - (g/(2*pi))*period^2*tanh(2*pi*depth/lambda);
end