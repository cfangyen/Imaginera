% % This script is for demontrating the process of using the software to
% analyze locomotor frequency from videos of moving C. elegans
clear;
% Locate the folder of the video data
ParentPath = uigetdir('D:\C elegans Team Dropbox\HONGFEI JI\Dian Hongfei shared\worm analysis\Data\High-Res Swim', 'Select parent folder');
% Initiation of the paramter for handling subfolders, files, and variables
sss = 1;
qqq = 1;
F_all = struct([]);
Fbout_all = struct([]);
groupDir = dir(fullfile(ParentPath, 'Exp_*'));
NumGroup = numel(groupDir);
% Set the framerate (fps) of the videos (metadata available in the property
% info of each video file)
fps = 30;
numPoints4analysis = round(fps*5);
NNumT = 0;
tic
for i = 1 : NumGroup % 
    currGroupName = fullfile(ParentPath, groupDir(i).name);
    expDir = dir(fullfile(currGroupName, '*_PCA_subtracted'));
    NumExp = numel(expDir);
    for j = 1 : NumExp
        currExpName = fullfile(currGroupName, expDir(j).name);
        cd(currExpName)
        Viddir = dir('*pca_*.avi');
        numVid = numel(Viddir);
        Fdata = struct([]);
        movmean_windowratio = 1/6;
        do_plot = 1;
        nbins = 16; % 16, 0.827
        delta = 0.045; % 0.045, 0.827
        ThresholdDelta = 0.8;
        MaxAreaVariation = 1;
        RegionAreaRange = [20 300];
        ttt = 1;
        for jj = 1 : numVid
            currVidname = Viddir(jj).name;
            fprintf('Analyzing Video: %s\n', currVidname)
            vid = VideoReader(currVidname);
            FrameRate = vid.FrameRate;
            AllFrames = read(vid);
            if contains(vid.VideoFormat, 'RGB')
                AllFrames = uint8(squeeze(mean(AllFrames, 3)));
            elseif contains(vid.VideoFormat, 'Grayscale')
                AllFrames = unit8(AllFrames);
            end
            NumT = size(AllFrames,3);
            AllFrames0 = AllFrames;
%             for kk = 1 : NumT
%                 currimg = AllFrames(:,:,kk);
%                 currimg = imresize(currimg, 0.1);
%                 AllFrames0 = cat(3, AllFrames0, currimg);
%             end
            NNumT = NNumT + NumT;
            if ~do_quickanalysis
                [F_avg, F_L, F_U, F_seq] = Frequency_analysis(AllFrames0, FrameRate, movmean_windowratio, do_plot, nbins,...
                                                              ThresholdDelta, MaxAreaVariation, RegionAreaRange, delta);
                Fdata(ttt).name  = currVidname;
                Fdata(ttt).F_seq = F_seq;
                Fdata(ttt).F_avg = F_avg;
                Fdata(ttt).F_L   = F_L;
                Fdata(ttt).F_U   = F_U;
                Fdata(ttt).F_gt  = FrameRate./median(diff(GT(ttt).BendTimes));
                Fdata(ttt).igt   = GT(ttt).BendTimes;
            else
                F_seq = Fdata(ttt).F_seq;
                F_gt  = Fdata(ttt).F_gt;
                igt   = Fdata(ttt).igt;
                [F_avg, F_L, F_U, imax, imin] = Quick_Frequency_analysis(F_seq, FrameRate, delta);
                nn = floor(numel(F_seq)/numPoints4analysis);
                for ii = 1 : nn
                    xL = 1+(ii-1)*numPoints4analysis;
                    if ii<nn
                        xR = 1+(ii+0)*numPoints4analysis;
                    else
                        xR = numel(F_seq);
                    end
                    imax_bout = imax(imax>=xL & imax<=xR);
                    imin_bout = imin(imin>=xL & imin<=xR);
                    igt_bout  = igt(igt>=xL & igt<=xR);
%                     if Fbout_avg/Fbout_gt<1.35
                        Fbout_all(qqq).F_avg = FrameRate./mean(cat(2, diff(imax_bout), diff(imin_bout)));
                        Fbout_all(qqq).F_gt  = FrameRate./mean(diff(igt_bout));
                        qqq = qqq+1;
%                     end
                end
            end
            F_all(sss).name  = currVidname;
            F_all(sss).F_seq = F_seq;
            F_all(sss).F_avg = F_avg;
            F_all(sss).F_L   = F_L;
            F_all(sss).F_U   = F_U;
            F_all(sss).F_gt  = FrameRate./mean(diff(GT(ttt).BendTimes));
%             fprintf('prd= %.2f Hz, gt= %.2f Hz\n', F_avg, FrameRate./mean(diff(GT(ttt).BendTimes)))
            ttt = ttt + 1;
            sss = sss + 1;
        end
        save('Fdata.mat', "Fdata")
    end
end
toc
cd(ParentPath)
save('Fdata_all_IJ_6bins_peakmedian.mat', 'F_all')

%% Plot Prediction vs GroundTruth (whole seq)
F_PRD = cat(2, F_all.F_avg);
F_GT  = cat(2, F_all.F_gt);
% jitterAmount = 0.05; % adjust this value as needed
% F_PRD_jt = F_PRD + jitterAmount * (rand(size(F_PRD)) - 0.5);
% F_GT_jt  = F_GT + jitterAmount * (rand(size(F_GT)) - 0.5);

CC = linspecer(4);
figure(1); clf
plot(F_GT, F_PRD, '.', 'Color',CC(1,:),'MarkerSize',6);
% scatter(F_GT, F_PRD, 15, 'filled')
hold on;
RMSE = sqrt(mean((F_GT - F_PRD).^2, 'omitnan'));
R2   = 1 - sum((F_GT - F_PRD).^2, 'omitnan')/sum((F_GT - mean(F_GT, 'omitnan')).^2, 'omitnan');
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

%% Plot Prediction vs GroundTruth (bout seq)
F_PRD = cat(2, Fbout_all.F_avg);
F_GT  = cat(2, Fbout_all.F_gt);
% % jitterAmount = 0.00; % adjust this value as needed
% F_PRD_jt = F_PRD + jitterAmount * (rand(size(F_PRD)) - 0.5);
% F_GT_jt  = F_GT + jitterAmount * (rand(size(F_GT)) - 0.5);

figure(3); clf
scatter(F_GT, F_PRD, 15, 'filled')
hold on;
RMSE = sqrt(mean((F_GT - F_PRD).^2));
R2   = 1 - sum((F_GT - F_PRD).^2)/sum((F_GT - mean(F_GT)).^2);
% Plot the line of perfect prediction
plot([min(F_GT), max(F_GT)], [min(F_GT), max(F_GT)], 'r', 'LineWidth', 2);
xlabel('Ground-truth undulatory frequency (Hz)');
ylabel('Predicted frequency (Hz)');
% Annotate the plot with the RMSE and R-squared values
text(2, 1.2, sprintf('RMSE: %.3f\nR^2: %.3f', RMSE, R2), 'FontSize', 16)
% Add a legend
legend('Data', 'x = y', 'Location', 'best');

hold off;

set(gca, 'FontSize', 16)

xlim([.8 3])
ylim([.8 3])