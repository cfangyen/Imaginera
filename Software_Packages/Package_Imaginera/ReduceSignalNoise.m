function V_filtered = ReduceSignalNoise(V, fps, movmean_windowratio)
% Reduce noise
timefilter = max([1 round(fps*movmean_windowratio)]);
filter1D = ones(1,timefilter)/timefilter;
V_median = medfilt1(V, timefilter, 'omitnan');
tmp = V_median;
V05 = prctile(tmp, 5);
V95 = prctile(tmp, 95);
V(V > V95) = V_median(V > V95);
V(V < V05) = V_median(V < V05);
V_padded = padarray(V, [0 timefilter-1], 'replicate', 'both');
V_filtpadded = conv(V_padded, filter1D, 'same');
V_filtered = V_filtpadded(timefilter:end-timefilter+1);


end

