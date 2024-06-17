%% Determine ROI of Wells
clc; clear; close all
NumWormsPerWell = 1;
ParentPath = uigetdir('D:\Hongfei-Dataset2024\Multiwell_Imaging\AD_expts\Treatment', 'Select parent folder');
originalVideoPath = 'original';
cd(ParentPath)
vidDir = dir(fullfile(originalVideoPath, '*.avi'));
NumVids = numel(vidDir);
sampleVidname = fullfile(vidDir(1).folder, vidDir(1).name);
sampleVid = VideoReader(sampleVidname);
ImgEx = im2gray(read(sampleVid, 1));
[ImgYsize, ImgXsize] = size(ImgEx);
figure(1); clf
set(gcf, 'Position', [72,162,1319,815])
h = imagesc(ImgEx); colormap gray; hold on;
set(gca,'XTick',[],'YTick',[],'FontSize', 20);
[~,~,FileType] = fileparts(sampleVidname);
Ans = inputdlg({sprintf('Selection Method\n1. Auto Selection (Web Method)\n2. Load ROI File')},'ROI Selection',1,{'1'});
SelectMethod = str2double(Ans{1});
switch SelectMethod
    case 1
        Fields = {'Number of Rows', 'Number of Columns'};
        Ans = inputdlg(Fields, 'ROI Auto-Selection', 1, {'8', '12'});
        NumRow = str2double(Ans{1});
        NumCol = str2double(Ans{2});
        NumROI = NumRow * NumCol;
        
        h.Parent.Title.String = 'Draw a circle to fit the UPPER-LEFT well.';
        hCorner = drawcircle('Color','r','FaceAlpha',0.4);
        Center = hCorner.Center;
        ULx = Center(1);
        ULy = Center(2);
        ULr = hCorner.Radius;

        h.Parent.Title.String = 'Draw a circle to fit the UPPER-RIGHT well.';
        hCorner = drawcircle('Color','r','FaceAlpha',0.4);
        Center = hCorner.Center;
        URx = Center(1);
        URy = Center(2);
        URr = hCorner.Radius;

        h.Parent.Title.String = 'Draw a circle to fit the LOWER-LEFT well.';
        hCorner = drawcircle('Color','r','FaceAlpha',0.4);
        Center = hCorner.Center;
        LLx = Center(1);
        LLy = Center(2);
        LLr = hCorner.Radius;

        h.Parent.Title.String = 'Draw a circle to fit the LOWER-RIGHT well.';
        hCorner = drawcircle('Color','r','FaceAlpha',0.4);
        Center = hCorner.Center;
        LRx = Center(1);
        LRy = Center(2);
        LRr = hCorner.Radius;

        ROI = struct([]);
        Lx = ULx + (LLx-ULx)/(NumRow-1)*(0:(NumRow-1));
        Ly = ULy + (LLy-ULy)/(NumRow-1)*(0:(NumRow-1));
        Rx = URx + (LRx-URx)/(NumRow-1)*(0:(NumRow-1));
        Ry = URy + (LRy-URy)/(NumRow-1)*(0:(NumRow-1));
        Radius = mean([ULr URr LLr LRr]);
        
        n = 0;
        for r = 1 : NumRow
            for c = 1 : NumCol
                n = n + 1;
                Center = [Lx(r)+(Rx(r)-Lx(r))/(NumCol-1)*(c-1) Ly(r)+(Ry(r)-Ly(r))/(NumCol-1)*(c-1)];
                ROI(n).Center = Center;
                ROI(n).Radius = Radius;
            end
        end
        
        for n = 1 : NumROI
            Center = ROI(n).Center;
            Radius = ROI(n).Radius;
            plot(Center(1),Center(2),'*r')
            drawcircle('Center', Center, 'Radius', Radius);
        end

        plot(ULx,ULy,'*m')
        plot(URx,URy,'*m')
        plot(LLx,LLy,'*m')
        plot(LRx,LRy,'*m')

        save(fullfile(ParentPath, 'ROI_Info.mat'),"ROI","NumRow","NumCol" )
    case 2
        ROI_file = uigetfile(ParentPath, 'Select ROI_Info.mat file');
        load(ROI_file, "ROI","NumRow","NumCol");
        NumROI = numel(ROI);
        for n = 1 : NumROI
            Center = ROI(n).Center;
            Radius = ROI(n).Radius;
            plot(Center(1),Center(2),'*r')
            drawcircle('Center', Center, 'Radius', Radius);
        end
end
fprintf('ROI loaded successfully!\n')
pause(1)
close

%% Initializing data analysis
happy = 0;
while happy==0
    fields={'Start Video',...
            'End Video',...
            'Analyze Frequency?(0: NO; 1: YES)',...
            'Frame/Time (sec)',...
            sprintf('Video File Prefix\nNote: All video files must be named \nin the format PREFIX####.ext'),...
            'BKG Threshold for ActVal Analysis',...
            'Spatial Filter Size',...
            'Frame Skip',...
            'Autosave per video?(0: NO; 1: YES)'};
    try
        Ans = inputdlg(fields,'Set Parameters',1,{num2str(NumStart),num2str(NumEnd),num2str(do_Frequency),...
                                                  num2str(fps),FilePrefix,num2str(NoiseThres),...
                                                  num2str(GaussianStd),num2str(frameSkip),num2str(AutoSave)});
    catch
        Ans = inputdlg(fields,'Set Parameters',1,{'1',num2str(NumVids),'1','5','','0.6','1.0','5','1'});
    end
    NumStart = str2double(Ans{1});
    NumEnd   = str2double(Ans{2});
    do_Frequency = str2double(Ans{3});
    fps = str2double(Ans{4});
    FilePrefix  = Ans{5};
    NoiseThres  = str2double(Ans{6});
    GaussianStd = str2double(Ans{7});
    frameSkip   = str2double(Ans{8});
    AutoSave    = str2double(Ans{9});
    
    % Calculate activity from two images to display for user
    x = -5:5;
    y = x;
    [xx, yy] = meshgrid(x, y);
    gau = exp(-sqrt(xx.^2+yy.^2)/GaussianStd^2);
    sample1 = double(im2gray(read(sampleVid, 1)));
    sample2 = double(im2gray(read(sampleVid, 1 + fps)));
    Diff1 = abs(sample1 - sample2)./(abs(sample1+sample2));
    smoothedDiff1 = convn(Diff1, gau, 'same');
    BinaryDiff1 = smoothedDiff1>NoiseThres;
    imagesc(Diff1, [0 0.7]); colormap('gray')
    set(gcf, 'Position', [72,162,1319,815])
    an = inputdlg('Click ok to display thresholded image');
    imagesc(BinaryDiff1); colormap('gray')
    myans=inputdlg({'Are you happy with parameters?  Type 1 if yes or 0 to reset: '},'Accept Parameters?',1,{'0'});
    happy=str2double(myans{1});
end

% AutoSave Path
if AutoSave ~= 0
    [SaveName, SavePath] = uiputfile([FilePrefix '_Result_DATA.mat'], 'AutoSave To');
else
    SaveName = [];
    SavePath = [];
end
fprintf('Parameter set successfully!\n')
pause(1)
close
NumVid = NumEnd-NumStart+1;
%% Output Variables Initialization
S = repmat(struct('Vidname', FilePrefix,...
                  'FrameRate', fps,...
                  'FrameSkip', frameSkip,...
                  'Freq', nan(NumRow, NumCol),...
                  'Freq_L', nan(NumRow, NumCol),...
                  'Freq_U', nan(NumRow, NumCol),...
                  'Fullseq', cell(NumRow, NumCol),...
                  'Periodicity', nan(NumRow, NumCol),...
                  'date', []), NumVid, 1);
NumFramesAll = zeros(NumVid, 1);

%% Main Frequency Analysis
if do_Frequency
    x = -5:5;
    y = x;
    [xx, yy] = meshgrid(x,y);
    gau = exp(-sqrt(xx.^2+yy.^2)/GaussianStd^2);
    [imxx,imyy] = meshgrid(1:ImgXsize,1:ImgYsize);
    issave_video = 0;

    movmean_windowratio = 1/6;
    do_plot = 0;
    nbins = 6;
    ThresholdDelta = 1.8;
    MaxAreaVariation = 0.4;
    RegionAreaRange = [20 300];
    % Make a folder for temporarily storing well video for compute frequency
    PathPCAWell = fullfile(ParentPath, 'PCA subtracted video of a well');
    if ~exist(PathPCAWell, 'dir')
        mkdir(PathPCAWell)
    end

    for i = 1:NumVids
        tic
        NumCurr = NumStart+i-1;
        NameCurr = [FilePrefix num2str(NumCurr, '%04d') FileType];
        CurrVid = VideoReader(fullfile(ParentPath, originalVideoPath, NameCurr));
%         NumFrames = 97;
        NumFrames = CurrVid.NumFrames;
        FrameRate = fps;
        % Preload Video Frames ---------------------------------------
        fprintf('Start preloading video: %s \n', NameCurr)
        AllFrames = read(CurrVid, [1 NumFrames]);
%         AllFrames = read(CurrVid);
        if contains(CurrVid.VideoFormat, 'RGB')
            AllFrames = squeeze(mean(AllFrames, 3));
        elseif contains(CurrVid.VideoFormat, 'Grayscale')
            AllFrames = squeeze(AllFrames);
        end
%         AllFrames = gpuArray(AllFrames);       
        fprintf('Complete preloading video: %s \n', NameCurr)
        ActVal = (S(i).ActVal)'/(400*2);
        ActVal = ActVal(:);
        Freq   = nan(NumROI, 1);
        Freq_L = nan(NumROI, 1);
        Freq_U = nan(NumROI, 1);
        Fullseq = cell(NumROI, 1);
        Periodicity = nan(NumROI, 1);
        WellVid_2d_All = cell(NumROI, 1);
        % PCA bkg subtraction and prestore all the well videos
%         parpool(3)
        for n = 1:NumROI
            % Determine ROI of the well
            Center = round(ROI(n).Center);
            xc = Center(1); yc = Center(2);
            R = round(ROI(n).Radius);
%             mask = false(ImgYsize, ImgXsize);
%             mask = mask | hypot(imxx-xc, imyy-yc) <= R;
%             mask_stacked = repmat(mask, [1 1 NumFrames]);
            H = 2*R+1;
            W = 2*R+1;
%             WellVid_2d = zeros(H*W, NumFrames);
            fprintf('PCA subtraction on well %d / %d of video %d / %d. \n', n, NumROI, i, NumVid)
%             AllFrames(~mask_stacked) = 0;
            WellVid_3d = AllFrames(yc-R:yc+R, xc-R:xc+R,:);
%             for jj = 1:NumFrames
%                 currimg = AllFrames(:,:,jj);
% %                 currimg(~mask) = 0;
%                 currimg = currimg(yc-R:yc+R, xc-R:xc+R);
%                 WellVid_2d(:, jj) = currimg(:);
%             end
            WellVid_2d = double(reshape(WellVid_3d, [H*W, NumFrames]));
            [score, first_pc] = pcasecon(WellVid_2d', 1);
%             [first_pc, score] = pca(WellVid_2d', 'Centered', false, 'NumComponents', 1);

            % The first principal component
            % first_pc = coeff(:,1);
            WellVid_2d_All{n} = WellVid_2d - first_pc' * score';
        end
        fprintf('Complete preloading all well videos: %s \n', NameCurr)
        clear AllFrames WellVid_2d WellVid_3d

        for n = 1:NumROI
            WellVid_2d = WellVid_2d_All{n};
            R = round(ROI(n).Radius);
            H = 2*R+1;
            W = 2*R+1;
            if issave_video
                % Write video
                vidname = [FilePrefix num2str(NumCurr, '%04d') sprintf('_Well%02d', n) FileType];
                WellVid_4d = reshape(WellVid_2d, [height, width, 1, NumFrames]);
                WellVid_4d = uint8(WellVid_4d);
                vidpath2save = fullfile(PathPCAWell, vidname);
                v = VideoWriter(vidpath2save);
                v.FrameRate = FrameRate;
                open(v);
                writeVideo(v, WellVid_4d)
                close(v)
                clear WellVid_4d v
            end
%             CurrWellVid = VideoReader(vidpath2save);
            fprintf('Analyzing well %d / %d of video %d / %d. \n', n, NumROI, i, NumVid)
            Period_score = ActVal(n);
            [WellF, WellF_L, WellF_U, WellFullseq]  = Frequency_analysis_mw(uint8(reshape(WellVid_2d, [H, W, NumFrames])),...
                                                          FrameRate, movmean_windowratio, do_plot, nbins,...
                                                          ThresholdDelta, MaxAreaVariation, RegionAreaRange);
%             clear WellVid_3d
%             fprintf('Frequency: %.2f Hz (%.2f~%.2f), ActLevel: %.3f\n', WellF, WellF_L, WellF_U, min([Period_score, 1]))
            Freq(n)   = WellF;
            Freq_L(n) = WellF_L;
            Freq_U(n) = WellF_U;
            Fullseq{n} = WellFullseq;
            Periodicity(n) = min([Period_score, 1]);
        end
        clear WellVid_2d_All
%         delete(gcp('nocreate'));
        toc
        S(i).Freq   = reshape(Freq,   [NumCol, NumRow])';
        S(i).Freq_L = reshape(Freq_L, [NumCol, NumRow])';
        S(i).Freq_U = reshape(Freq_U, [NumCol, NumRow])';
        S(i).Fullseq = reshape(Fullseq, [NumCol, NumRow])';
        S(i).Periodicity = reshape(Periodicity, [NumCol, NumRow])';
        save(fullfile(SavePath,SaveName), 'S')
        if issave_video
            % Clear the PCA subtracted video folder to perform frequency analysis of the next video
            all_items = dir(PathPCAWell); % List all files and directories
            % Loop through the list
            for k = 1:length(all_items)
                item_name = all_items(k).name;
                % Skip the current folder '.' and parent folder '..'
                if strcmp(item_name, '.') || strcmp(item_name, '..')
                    continue;
                end
                item_path = fullfile(PathPCAWell, item_name);
                if all_items(k).isdir  % If it's a directory
                    rmdir(item_path, 's'); % Remove directory and its contents
                else
                    delete(item_path); % Remove file
                end
            end
        end
    end
end
