function [F_avg, F_L, F_U, F_seq]...
         = Frequency_analysis(Vid3d, fps, movmean_windowratio, do_plot, nbins, ThresholdDelta, MaxAreaVariation, RegionAreaRange, delta)
%% Hu analysis of video 
SScale = 1; TScale = 1;
I = Hu_calculation(Vid3d, SScale, TScale);
[~, I_nr] = Hu_feature_extraction(I, fps, movmean_windowratio);
%% MSER analysis of video
J_feature = MSER_calculation(Vid3d, SScale, TScale,...
                             ThresholdDelta, MaxAreaVariation, RegionAreaRange, do_plot);
[~, J_nr] = MSER_feature_extraction(J_feature, fps, nbins, movmean_windowratio);
%% Combine Hu and SIFT analysis
IJ_nr = cat(2, I_nr, J_nr);
%% Frequency analysis
AA = cosineSimilarity(IJ_nr);
cov_IJ_nr = AA;
[F_avg, F_L, F_U, F_seq] = Frequency_extraction1(cov_IJ_nr, fps, delta, movmean_windowratio);
end

