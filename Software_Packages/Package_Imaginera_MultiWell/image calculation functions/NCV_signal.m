function NCV = NCV_signal(signal)
%NCV_SIGNAL Measure the variability of a time sequence
% Calculate mean and standard deviation
n = numel(signal);
mu = mean(signal, 'omitnan');
sigma = std(signal, 'omitnan');

% Coefficient of Variation
CV = sigma / mu;

% Maximum possible standard deviation for a binary signal with mean 0.5
signal_max_var = [zeros(round(n/2), 1); ones(round(n/2), 1)]; % half zeros, half ones
sigma_max = std(signal_max_var);
CV_max = sigma_max / 0.5;

% Normalized Coefficient of Variation
NCV = CV / CV_max;
end

