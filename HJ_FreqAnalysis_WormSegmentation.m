% % This script is for demontrating the process of using the segmentation 
% software to analyze locomotor frequency from videos of moving C. elegans
clc; close all; clear
% Locate the folder of the video data
TargetFolder = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol\Noisy videos\p0243';

ParentPath = uigetdir(TargetFolder, 'Select parent folder');
% global start_illum end_illum prefix pathname filename
% global conc pix_per_mm wormthreshold isie decim filsize spline_p initials
qqq = 1;
cd(ParentPath)
load('GTdata_curv.mat', 'GTdata')
load('CVdata_curv.mat', 'CVdata')
viddir = dir('*.avi');
NumVid = numel(viddir);
SScale = [0.03, 0.1, 0.3, 1];
TScale = [1/15 1/10 1/6 1/3 1/2 1];
% Set the framerate (fps) of the videos (metadata available in the property
% info of each video file)
do_quickanalysis = 0;
sample_vid = VideoReader(viddir(1).name);
fps = sample_vid.FrameRate;
numPoints4analysis = round(fps*5);
NNumT = 0;
movmean_windowratio = 1/6;
do_plot = 0;
nbins = 16; % 16, 0.827

resizefactor = SScale(end);
tscale = TScale(end);
% Parameter initialization

conc = 0;
wormlabel = 1;
pix_per_mm = 1;
wormthreshold = 0.10;
decim = 1;
filsize = 0.2;
start_illum = 1;
end_illum = 1;
spline_p = 0.01;
delta = 0.05;

seg_info = struct([]);

for i = 1 : NumVid
    currVidname = viddir(i).name;
    fprintf('Analyzing Video: %s\n', currVidname)
    vidObj = VideoReader(currVidname);
    NumFrames = vidObj.NumFrames;
    fps       = vidObj.FrameRate;

    isie = [1, NumFrames];

    %  Compute undulatory variables
    do_const = 0;

    options = {conc, wormlabel, fps, pix_per_mm, wormthreshold,...
        isie, decim, filsize, start_illum, end_illum, spline_p};

    curv = WORMSHAPE_MAINCALCULATION(vidObj, options, resizefactor, tscale);
    figure(1); clf
    imagesc(curv);
    if ~do_quickanalysis
        [F_avg, F_L, F_U, F_seq] = Frequency_extraction_curv(curv, fps, delta);
    else
    end
    CVdata(i).name  = currVidname;
    CVdata(i).F_seq = F_seq;
    CVdata(i).F_avg = F_avg;
    CVdata(i).F_L   = F_L;
    CVdata(i).F_U   = F_U;

    Fgt = GTdata(i).F_gt_mean;
    Fpr = F_avg;
    pct_change = abs(Fpr-Fgt)/Fgt;
    if pct_change>=.1
        CVdata(i).isoutlier = 1;
    else
        CVdata(i).isoutlier = 0;
    end
    fprintf('Predict: %.3f Hz, GrTruth: %.3f Hz, isoutlier: %d\n', Fpr, Fgt, CVdata(i).isoutlier)
    fprintf('\n')
    CVdata(i).F_gt  = Fgt;
end

F_PRD = cat(2, CVdata.F_avg);
F_GT  = cat(2, CVdata.F_gt);
% jitterAmount = 0.05; % adjust this value as needed
% F_PRD_jt = F_PRD + jitterAmount * (rand(size(F_PRD)) - 0.5);
% F_GT_jt  = F_GT + jitterAmount * (rand(size(F_GT)) - 0.5);

CC = linspecer(4);
figure(1); clf
plot(F_GT, F_PRD, '.', 'Color',CC(1,:),'MarkerSize',6);
% scatter(F_GT, F_PRD, 15, 'filled')
hold on;
RMSE = sqrt(mean((F_GT - F_PRD).^2, 'omitnan'));
disp(['RMSE: ', num2str(RMSE)]);
R2   = 1 - sum((F_GT - F_PRD).^2, 'omitnan')/sum((F_GT - mean(F_GT, 'omitnan')).^2, 'omitnan');
disp(['Regular R²: ', num2str(R2)]);
mdl_robust = fitlm(F_GT, F_PRD, 'RobustOpts', 'on');
R2_robust = mdl_robust.Rsquared.Ordinary;
disp(['Robust R²: ', num2str(R2_robust)]);
NAN  = sum(isnan(F_PRD));
disp(['# of datapoint excluded: ', num2str(NAN)]);
% Plot the line of perfect prediction
plot([min(F_GT), max(F_GT)], [min(F_GT), max(F_GT)], 'r', 'LineWidth', 2);
xlabel('Ground-truth (Hz)');
ylabel('Automated (Hz)');
% Annotate the plot with the RMSE and R-squared values
% text(2, 1.2, sprintf('RMSE: %.3f\nR^2: %.3f', RMSE, R2), 'FontSize', 16)
% Add a legend
% legend('Data', 'x = y', 'Location', 'best');

hold off;

set(gca, 'FontSize', 6, 'FontName', 'Arial', 'Box','off','DataAspectRatio',[1 1 1])
xticks([1 1.5 2 2.5])
yticks([1 1.5 2 2.5])
% xticklabels({'1.0', '1.5', '2.0', '2.5'})
yticklabels({'1.0', '1.5', '2.0', '2.5'})
xlim([1 2.7]);
ylim([1 2.7]);

set(gcf, 'Position', [295,754,110,110])
save('CVdata_curv.mat', 'CVdata')