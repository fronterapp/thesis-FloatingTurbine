%% Computes phase averaged value of a periodic signal, refFun
%  Computation will not take into account the first N = delayPeriods periods
%  Each period will be shortened or enlarged up to M = refSampleLength points

function [phaseAvg] = averagePhase(ref_fun, delay_periods, ref_sample_length)

    % Find the peaks
    [peaks, idx]= findpeaks(ref_fun);
    % To avoid numerical error, consider peaks only those that are at least 5% of wave height
    idxPeaks = idx(peaks>0.05*max(peaks));
    % Split according to the peaks' position
    splitFun = zeros(length(idxPeaks)-1-delay_periods,ref_sample_length); 

    for i=delay_periods+1:length(idxPeaks)-1
        % Size of the actual period
        sampleLength = idxPeaks(i+1)-idxPeaks(i)+1; 
        % Interpolate in case the periods do not contain the same number of samples, needed if we want to sum them up
        splitFun(i,:) = interp1(linspace(0,1,sampleLength),ref_fun(idxPeaks(i):idxPeaks(i+1)),linspace(0,1,ref_sample_length));
    end

    % Compute the phase-average
    phaseAvg = sum(splitFun,1)./(length(idxPeaks)-1-delay_periods);

end