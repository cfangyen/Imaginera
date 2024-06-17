function [F_avg, F_L, F_U, Fullseq]...
         = Frequency_analysis_mw(WellVid_3d, fps, movmean_windowratio, do_plot, nbins,ThresholdDelta, MaxAreaVariation, RegionAreaRange)
%% Hu analysis of video 
SScale = 1; TScale = 1;
I = Hu_calculation_mw(gpuArray(WellVid_3d), SScale, TScale);
[~, I_nr] = Hu_feature_extraction_mw(I, fps, movmean_windowratio);
timefilter = max(1, round(fps*movmean_windowratio));
covfilter = fspecial('average',[timefilter,timefilter]);

%% MSER analysis of video
J_feature = MSER_calculation_mw(WellVid_3d, SScale, TScale,...
                                ThresholdDelta, MaxAreaVariation, RegionAreaRange, do_plot);
[~, J_nr] = MSER_feature_extraction_mw(J_feature, fps, nbins, movmean_windowratio);
%% Combine Hu and SIFT analysis
IJ_nr = cat(2, I_nr, J_nr);
%% Frequency analysis
% cov_I_nr = cov(I_nr', 'includenan');
% cov_I_filtered = imfilter(cov_I_nr, covfilter, 'replicate');
% cov_J_nr = cov(J_nr', 'includenan');
% cov_J_filtered = imfilter(cov_J_nr, covfilter, 'replicate');
cov_IJ_nr = cov(IJ_nr', 'partialrows');
cov_IJ_filtered = imfilter(cov_IJ_nr, covfilter, 'replicate');
[F_avg, F_L, F_U, Fullseq] = Frequency_extraction_mw(cov_IJ_filtered, fps);
%%
% NumT = size(IJ_nr, 1); T = (0:(NumT-1))/5;
% m=100;
% % cm_parula  =fake_parula(m);
% % cm_magma   =magma(m);
% cm_inferno =inferno(m);
% % cm_plasma  =plasma(m);
% % cm_viridis =viridis(m);
% figure(5); clf
% imagesc(T,T,cov_IJ_filtered);
% colormap(cm_inferno)
% % clim([0 1])
% xticks([0 T(end)])
% yticks([0 T(end)])
% ylabel('Time (s)')
% xlabel('Time (s)')
% set(gca, 'FontSize', 6, 'FontName','Arial','PlotBoxAspectRatio', [1 1 1])
% set(gcf, 'Position', [295,754,95,95])
% 
% figure(11); clf
% plot(T, Fullseq)
% xlabel('Time (s)')
% drawnow

end

