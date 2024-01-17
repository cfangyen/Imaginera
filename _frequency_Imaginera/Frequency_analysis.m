function [F_avg, F_L, F_U, F_seq]...
         = Frequency_analysis(Vid3d, fps, movmean_windowratio, do_plot, nbins, ThresholdDelta, MaxAreaVariation, RegionAreaRange, delta)
%% Hu analysis of video 
SScale = 1; TScale = 1;
I = Hu_calculation(Vid3d, SScale, TScale);
[~, I_nr] = Hu_feature_extraction(I, fps, movmean_windowratio);
timefilter = max(1, round(fps*movmean_windowratio));
covfilter = fspecial('average',[timefilter,timefilter]);

%% MSER analysis of video
J_feature = MSER_calculation(Vid3d, SScale, TScale,...
                             ThresholdDelta, MaxAreaVariation, RegionAreaRange, do_plot);
[~, J_nr] = MSER_feature_extraction(J_feature, fps, nbins, movmean_windowratio);
%% Combine Hu and SIFT analysis
IJ_nr = cat(2, I_nr, J_nr);
%% Frequency analysis
% cov_J_nr  = cov(J_nr', 'includenan');
% cov_IJ_nr = cov(IJ_nr', 'includenan');
% cov_IJ_nr = xcorr(IJ_nr', 'normalized');
% AA = cov_IJ_nr(size(IJ_nr,2),:);
% AA = reshape(AA, [size(IJ_nr,1) size(IJ_nr,1)]);
IJ_nr = IJ_nr(510:end, :);
IJ_nr = zscore(IJ_nr,0, 1);
AA = cosineSimilarity(IJ_nr);
cov_IJ_nr = AA;
N = size(cov_IJ_nr, 1);
auto_IJ = zeros(N, 1);
for i = 1 : N
    auto_IJ(i) = mean(diag(cov_IJ_nr, i-1), 'omitnan');
end
figure(1); clf
plot((0: numel(auto_IJ)-1)/30, auto_IJ)
xlim([0 3])
xlabel('Time (s)')
ylabel("Autocorrelation of K")
set(gca, 'FontSize', 6, 'FontName', 'Arial', 'Box','off')

set(gcf, 'Position', [295,754,200,160])
% cov_J_filtered = imfilter(cov_J_nr, covfilter, 'replicate');
% cov_IJ_filtered = imfilter(cov_IJ_nr, covfilter, 'replicate');
%%%%%%%%%%%%
% AA = cov_IJ_nr(size(IJ_nr,2),:);
% AA = reshape(AA, [size(IJ_nr,1) size(IJ_nr,1)]);
% % AA = imfilter(AA, covfilter, 'replicate');
% figure(1)
% imagesc(AA)
%%%%%%%%%%%%
[F_avg, F_L, F_U, F_seq] = Frequency_extraction(cov_IJ_nr, fps, delta);

%%
% %% large-amplitude movement (1-177)
% sumimg = zeros(size(Vid3d(:,:,1)));
% idx = 1:2:120;
% for i = 1:numel(idx)
%     currimg = double(Vid3d(:,:,idx(i))).*i^2;
%     sumimg = sumimg+currimg;
% end
% figure(6);clf
% imagesc(sumimg); axis image;
% colormap('turbo')
% set(gca, 'FontName', 'Arial', 'FontSize', 6)
% set(gcf, 'Position', [84,747,280,225])
% %% Constrained movement (178-456)
% sumimg = zeros(size(Vid3d(:,:,1)));
% idx = 178:2:456;
% for i = 1:numel(idx)
%     currimg = double(Vid3d(:,:,idx(i))).*i^2;
%     sumimg = sumimg+currimg;
% end
% figure(7);clf
% imagesc(sumimg); axis image;
% colormap('turbo')
% set(gca, 'FontName', 'Arial', 'FontSize', 6)
% set(gcf, 'Position', [84,747,280,225])
% %% Reverse and pirouette (456-621)
% sumimg = zeros(size(Vid3d(:,:,1)));
% idx = 456:2:570;
% for i = 1:numel(idx)
%     currimg = double(Vid3d(:,:,idx(i))).*i^2;
%     sumimg = sumimg+currimg;
% end
% figure(8);clf
% imagesc(sumimg); axis image;
% colormap('turbo')
% set(gca, 'FontName', 'Arial', 'FontSize', 6)
% set(gcf, 'Position', [84,747,280,225])
% %% Undulation while turning (622-834)
% sumimg = zeros(size(Vid3d(:,:,1)));
% idx = 622:2:834;
% for i = 1:numel(idx)
%     currimg = double(Vid3d(:,:,idx(i))).*i^2;
%     sumimg = sumimg+currimg;
% end
% figure(9);clf
% imagesc(sumimg); axis image;
% colormap('turbo')
% set(gca, 'FontName', 'Arial', 'FontSize', 6)
% set(gcf, 'Position', [84,747,280,225])
% %% Undulation at a steady orientation (835-1001)
% sumimg = zeros(size(Vid3d(:,:,1)));
% idx = 860:2:1001;
% for i = 1:numel(idx)
%     currimg = double(Vid3d(:,:,idx(i))).*i^2;
%     sumimg = sumimg+currimg;
% end
% figure(10);clf
% imagesc(sumimg); axis image;
% colormap('turbo')
% set(gca, 'FontName', 'Arial', 'FontSize', 6)
% set(gcf, 'Position', [84,747,280,225])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 100;
cm_magma=magma(m);
NumT = size(IJ_nr, 1); T = (0:(NumT-1))/30;
NumS = size(IJ_nr, 2); S = 1:NumS;
IJ_nr_normal = zeros(size(IJ_nr));
for i = 1 : NumS
    V  = IJ_nr(:, i);
    mu = mean(V);
    sd = std(V);
    Vn = (V-mu)/sd;
    IJ_nr_normal(:,i) = Vn;
end
figure(4); clf
imagesc(T,S,IJ_nr_normal')
hold on;
plot([0 T(end)], [8.5 8.5], '-w')
plot([0 T(end)], [24.5 24.5], '-w')
hold off;
colormap(cm_magma)
% colormap(cmap_redblue(0.7))
xlabel('Time (s)')
ylabel(sprintf('Invariant index'))
yticks([4.5 16.5 32.5])
yticklabels({'H','E','O'})
    clim([-2 4]);
%     xlim([15 25])
%     xlim([17.5 19])
%     colorbar
set(gcf, 'Position', [295,754,380,160])
set(gca, 'FontName', 'Arial', 'FontSize', 6)
drawnow

m=100;
cm_parula  =fake_parula(m);
cm_magma   =magma(m);
cm_inferno =inferno(m);
cm_plasma  =plasma(m);
cm_viridis =viridis(m);
figure(5); clf
imagesc(T,T,cov_IJ_nr);
colormap(cm_inferno)
% xticks([40 80 120 160])
% yticks([40 80 120 160])
ylabel('Time (s)')
xlabel('Time (s)')
% xlim([15 25])
% ylim([15 25])
set(gca, 'FontSize', 6, 'FontName','Arial','PlotBoxAspectRatio', [1 1 1])
set(gcf, 'Position', [295,754,240,240])

figure(11); clf
plot(T, -mean(cov_IJ_nr, 1, 'omitnan'))
% xlim([15 25])
xlabel('Time (s)')
drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual traces of feature invariants
ImagineraFeatures = IJ_nr(:, [1:8 (8+[1 6 11 16]) (24+[1 6 11 16])]);
% ImagineraFeatures = log(zscore(ImagineraFeatures, 0, 1)+10);
ImagineraFeatures(:,[1:8 13:16]) = zscore(ImagineraFeatures(:,[1:8 13:16]), 0, 1);
figure(12); clf;
hold on
plot(T, ImagineraFeatures(:,1))
plot(T, ImagineraFeatures(:,2))
plot(T, ImagineraFeatures(:,3))
plot(T, ImagineraFeatures(:,4))
xlim([17.5 19])
set(gcf, 'Position', [295,754,380,80])
set(gca, 'FontSize', 6, 'FontName','Arial')

figure(13); clf;
hold on
plot(T, ImagineraFeatures(:,5))
plot(T, ImagineraFeatures(:,6))
plot(T, ImagineraFeatures(:,7))
plot(T, ImagineraFeatures(:,8))
xlim([17.5 19])
set(gcf, 'Position', [295,754,380,80])
set(gca, 'FontSize', 6, 'FontName','Arial')

figure(14); clf
hold on
plot(T, zscore(log(gradient(ImagineraFeatures(:,9))+1)))
plot(T, zscore(log(gradient(ImagineraFeatures(:,10))+1)))
plot(T, zscore(log(gradient(ImagineraFeatures(:,11))+1)))
plot(T, zscore(log(gradient(ImagineraFeatures(:,12))+1)))
xlim([17.5 19])
set(gcf, 'Position', [295,754,380,80])
set(gca, 'FontSize', 6, 'FontName','Arial')

% figure(14); clf
% hold on
% plot(T, -(log(-sin(ImagineraFeatures(:,9))+1.1)+.7))
% plot(T, -(log(-sin(ImagineraFeatures(:,10))+1.1)+.7))
% plot(T, -(log(-sin(ImagineraFeatures(:,11))+1.1)+.7))
% plot(T, -(log(-sin(ImagineraFeatures(:,12))+1.1)+.7))
% set(gcf, 'Position', [295,754,380,80])

figure(15); clf
hold on
plot(T, ImagineraFeatures(:,13))
plot(T, ImagineraFeatures(:,14))
plot(T, ImagineraFeatures(:,15))
plot(T, ImagineraFeatures(:,16))
xlim([17.5 19])
set(gcf, 'Position', [295,754,380,80])
set(gca, 'FontSize', 6, 'FontName','Arial')

%%
end

