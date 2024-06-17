function [I, I_noiseremoved] = Hu_feature_extraction_mw(I, fps, movmean_windowratio)
%HU_FEATURE_EXTRACTION: extracting features from the raw Hu's invariants
I_noiseremoved = nan(size(I));
numHu = size(I, 2);
for i = 1 : numHu
    V = I(:, i);
    V_filtered = ReduceSignalNoise_mw(V', fps, movmean_windowratio);
    V_filtered = detrend(V_filtered);
    I_noiseremoved(:, i) = V_filtered;
end
end

