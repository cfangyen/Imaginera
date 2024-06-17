% % This script is for demontrating the process of using the Imaginera
% software to analyze locomotor frequency from videos of moving C. elegans
clear; clc; close all
%% Locate the folder of the video data
TargetFolder = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol\Noisy videos\p0000';

ParentPath = uigetdir(TargetFolder, 'Select parent folder');
%% Initializing data analysis
fields = {'Ground-truth data available?(0: NO; 1: YES)',...
          'Win-size for data smoothing(default:1/6 sec)',...
          'Do plot during analysis?(0: NO; 1: YES)',...
          '# of MSER regions(default:16)',...
          'MSER ThresholdDelta(default:1;higher value for strict worm detection)',...
          'MSER MaxAreaVariation(default:1;lower value for strict worm detection)',...
          'RegionAreaRange(default:[20 300];range should cover worm size)',...
          'Prominence for peak/trough detection(default:0.14)'};
try
    Ans = inputdlg(fields,'Set Parameters',1,{num2str(do_loadGT),num2str(movmean_windowratio),num2str(do_plot),...
                                               num2str(nbins),num2str(ThresholdDelta),num2str(MaxAreaVariation),...
                                               mat2str(RegionAreaRange),num2str(Prom)});
catch
    Ans = inputdlg(fields,'Set Parameters',1,{'1','1/6','0','16','1','1','[20 300]','0.14'});
end
% Locate the dataset folder
cd(ParentPath)
% Set parameters for data analysis
do_loadGT = str2double(Ans{1});
movmean_windowratio = str2double(Ans{2});
do_plot = str2double(Ans{3});
nbins = str2double(Ans{4});
ThresholdDelta = str2double(Ans{5});
MaxAreaVariation = str2double(Ans{6});
RegionAreaRange = Str2Mat(Ans{7});
Prom = str2double(Ans{8});
% Loading the ground-truth dataset
if do_loadGT
    load('GTdata.mat', 'GTdata')
end
load('CVdata.mat', 'CVdata')
viddir = dir('*.avi');
NumVid = numel(viddir);
tic
% Perform frequency measurement video by video
for i = 1 : NumVid
    % Fetch video name
    currVidname = viddir(i).name;
    fprintf('Analyzing Video: %s\n', currVidname)
    % Read video data
    vid = VideoReader(currVidname);
    % Acquire video framerate
    FrameRate = vid.FrameRate;
    % Assign video data to a matrix variable MxNxP
    AllFrames = read(vid);
    % Regularize data into grayscale format
    if contains(vid.VideoFormat, 'RGB')
        AllFrames = uint8(squeeze(mean(AllFrames, 3)));
    elseif contains(vid.VideoFormat, 'Grayscale')
        AllFrames = unit8(AllFrames);
    end
    % Number of frames of the video
    NumT = size(AllFrames,3);
    AllFrames0 = AllFrames;
    % Feed the video data into Imaginera algorithm to measure frequency
    [F_avg, F_L, F_U, F_seq] = Frequency_analysis(AllFrames0, FrameRate, movmean_windowratio, do_plot, nbins,...
        ThresholdDelta, MaxAreaVariation, RegionAreaRange, Prom(1));
    % Save video name
    CVdata(i).name  = currVidname;
    % Save the autocorrelation sequency
    CVdata(i).F_seq = F_seq;
    % Save the prediction of moving frequency
    CVdata(i).F_avg = F_avg;
    % Save the lower bound of frequency
    CVdata(i).F_L   = F_L;
    % Save the upper bound of frequency
    CVdata(i).F_U   = F_U;
    % Fetch the ground-truth of frequency for the current video
    if do_loadGT
        Fgt = GTdata(i).F_gt_mean;
    end
    Fpr = F_avg;
    % Comparison between prediction and ground-truth for the current video
    if do_loadGT
        pct_change = abs(Fpr-Fgt)/Fgt;
        if pct_change>=.1
            CVdata(i).isoutlier = 1;
        else
            CVdata(i).isoutlier = 0;
        end
        fprintf('Predict: %.3f Hz, GrTruth: %.3f Hz, isoutlier: %d\n', Fpr, Fgt, CVdata(i).isoutlier)
        CVdata(i).F_gt  = Fgt;
    else
        fprintf('Predict: %.3f Hz\n', Fpr)
    end
    fprintf('\n')
end

toc
% Save all the predictions into a dataset
save('CVdata.mat', 'CVdata')

%% Plot Prediction vs GroundTruth (whole seq)
if do_loadGT
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
end
