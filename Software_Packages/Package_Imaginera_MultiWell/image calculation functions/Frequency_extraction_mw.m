function [F0, F0_L, F0_U, Fullseq] = Frequency_extraction_mw(C, fps)
%FREQUENCY_EXTRACTION
% % Compute frequency using peakfinding after avg in 1 direction
delta = 0.04;
avgC2analyze = mean(C, 1, 'omitnan');
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
    %     [imax, imin] = C2_get_curvature_peaks(stdC2analyze, 0);
    %     [imax, ~] = verify_extrema(stdC2analyze, imax, imin);
    %     T_com = cat(1, diff(imax), diff(imin));
    %     T = median(T_com)/fps;
    %     T_L = prctile(T_com, 25)/fps;
    %     T_U = prctile(T_com, 75)/fps;
    %     F = 1/T;
    %     F_L = 1/T_U;
    %     F_U = 1/T_L;
else
    F0   = NaN;
    F0_L = NaN;
    F0_U = NaN;
end
Fullseq = avgC2analyze;

% if sum(isfinite(avgC2analyze))>=6
%     figure(2); clf
%     subplot(211)
%     findpeaks(avgC2analyze,  'MinPeakProminence', Prominence);
%     ylabel('AVG of COV1d')
%     subplot(212)
%     findpeaks(-avgC2analyze, 'MinPeakProminence', Prominence);
%     ylabel('AVG of COV1d')
%     xlabel('Frame index')
%     set(gca, 'FontSize', 12)
% end
% 
% figure(3); clf
% imagesc(C); colorbar;
% title(sprintf('F(Hz): %.2f, Fl: %.2f, Fu: %.2f', F0, F0_L, F0_U))
% ylabel('Frame index')
% xlabel('Frame index')
% set(gca, 'FontSize', 12)
% 
% drawnow
% 
% pause(0.001)

end

