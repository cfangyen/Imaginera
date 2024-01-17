function [F0, F0_L, F0_U, F_seq] = Frequency_extraction(C, fps, delta)
%FREQUENCY_EXTRACTION
% % Compute frequency using peakfinding after avg in 1 direction
% avgC2analyze = mean(C, 1, 'omitnan');
% avgC2analyze = C(1,:);
N = size(C, 1);
auto_IJ = zeros(N, 1);
for i = 1 : N
    auto_IJ(i) = mean(diag(cov_IJ_nr, i-1), 'omitnan');
end
avgC2analyze = auto_IJ;
% [VV, ~] = eig(C);
% VV(:, end)
% avgC2analyze = VV(:, end)';

if sum(isfinite(avgC2analyze))>=6
    avgmedian  = median(avgC2analyze-min(avgC2analyze), 'omitnan');
    Prominence = avgmedian*delta;
    [~, imax0] = findpeaks(avgC2analyze,  'MinPeakProminence', Prominence);
    [~, imin0] = findpeaks(-avgC2analyze, 'MinPeakProminence', Prominence);

    T0_com = cat(2, diff(imax0), diff(imin0));
    if numel(T0_com)>=2
        T0 = mean(T0_com)/fps;
        T0_L = prctile(T0_com, 25)/fps;
        T0_U = prctile(T0_com, 75)/fps;
        F0 = 1/T0;
        F0_L = 1/T0_U;
        F0_U = 1/T0_L;
    else
        F0 = NaN;
        F0_L = NaN;
        F0_U = NaN;
    end
else
    F0   = NaN;
    F0_L = NaN;
    F0_U = NaN;
end
F_seq = avgC2analyze;
% m=100;
% cm_viridis=viridis(m);
% figure(2); clf
% subplot(211)
% findpeaks(avgC2analyze,  'MinPeakProminence', Prominence);
% ylabel('AVG of COV1d')
% subplot(212)
% findpeaks(-avgC2analyze, 'MinPeakProminence', Prominence);
% ylabel('AVG of COV1d')
% xlabel('Frame index')
% set(gca, 'FontSize', 12)
% 
% figure(3); clf
% imagesc(C);
% colormap(cm_viridis)
% title(sprintf('F(Hz): %.2f, Fl: %.2f, Fu: %.2f', F0, F0_L, F0_U))
% % xticks([40 80 120 160])
% % yticks([40 80 120 160])
% ylabel('Frame index')
% xlabel('Frame index')
% set(gca, 'FontSize', 6, 'FontName','Arial','PlotBoxAspectRatio', [1 1 1])
% % set(gcf, 'Position', [295,754,180,180])
% 
% drawnow
% 
% pause(1)

end


