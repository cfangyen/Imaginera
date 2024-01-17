function Vid_BkgSubtraction_Output(filename, orgviddir, wellimgdir, wellviddir, fieldmode, numwells, separation)
%%
% Load video
parentFolder = fileparts(orgviddir);
[~, vidname] = fileparts(filename);
Vidfullname = fullfile(orgviddir, filename);
vidObj = VideoReader(fullfile(orgviddir, filename));
numFrames = ceil(vidObj.Duration * vidObj.FrameRate);
FrameSample = rgb2gray(read(vidObj,1));
FPS = vidObj.FrameRate;
[numRow, numCol] = size(FrameSample);
%% Draw ROIs for wells
k = 1;
% Identify wells
ff = figure('Position',[500 500 400 300]);
imshow(FrameSample)
ff.WindowState = 'maximized';
axis('on', 'image');
title(sprintf('Indicate %d wells', numwells), 'FontSize', 15);
maskfile = dir(fullfile(parentFolder, sprintf('masks%d.mat', numwells)));
if isempty(maskfile)
    roi = cell(numwells, 1);
    for kk = 1 : numwells
        uiwait(helpdlg('Please click and drag out an ROI.'));
        roi{kk} = drawcircle('Color','r','FaceAlpha',0.4);
        title('Double click')
        waitfordoubleclick
        disp('Program execution has resumed');
%         savepathname_vid = fullfile(wellviddir,sprintf('Well%d',kk));
        savepathname_img = fullfile(wellimgdir,sprintf('Well%d',kk));
%         mkdir(savepathname_vid);
        mkdir(savepathname_img);
    end
    save(fullfile(parentFolder, sprintf('masks%d.mat', numwells)),'roi');
    load(fullfile(parentFolder, sprintf('masks%d.mat', numwells)),'roi')
else
    for kk = 1 : numwells
%         savepathname_vid = fullfile(wellviddir,sprintf('Well%d',kk));
        savepathname_img = fullfile(wellimgdir,sprintf('Well%d',kk));
%         mkdir(savepathname_vid);
        mkdir(savepathname_img);
    end
    load(fullfile(parentFolder, sprintf('masks%d.mat', numwells)),'roi')
end
close;

%% Extract well images from videos
% Generate the background
fprintf('-------------------------------------Extract well images from videos-------------------------------------\n')
while k <= numFrames
    if ~hasFrame(vidObj)
        break
    end
    currentFrame = rgb2gray(readFrame(vidObj));
    if mod(k, FPS) == 0
        for ii = 1 : numwells
            angles = linspace(0, 2*pi, 10000);
            x = cos(angles) * roi{ii}.Radius + roi{ii}.Center(1);
            y = sin(angles) * roi{ii}.Radius + roi{ii}.Center(2);
            mask = poly2mask(x, y, numRow, numCol);
            maskedImage = currentFrame; % Initialize with the entire image.
            maskedImage(~mask) = nan; % Zero image outside the circle mask.
            props = regionprops(mask, 'BoundingBox');
            maskedImage = imcrop(maskedImage, props.BoundingBox);
            savepathname_img = fullfile(wellimgdir,sprintf('Well%d',ii));
            imgname = fullfile(savepathname_img, [vidname sprintf('_frame%d_W%02d.png', k, ii)]);
            imwrite(maskedImage, imgname);
        end
    end
    currentFrame = double(currentFrame);
    if k == 1
        background = currentFrame;
    else
        if fieldmode == 'd' % if dark-field
            background = min(background, currentFrame);
        elseif fieldmode == 'b' % if bright-field
            background = max(background, currentFrame);
        end
    end
    k = k + 1;
end
fprintf('-------------------------------------Finished-------------------------------------\n')
pause(2)
%% Calculate and store worm activity for each one-well videos
background = uint8(background);
fprintf('-------------------------Calculate and store worm activity-------------------------\n')
% delete(gcp('nocreate'))
% parpool(12)
[LocoAVG, LocoSEM] = Vid_locoactivity(Vidfullname, separation, background, roi);
% delete(gcp('nocreate'))
savename_matfile = fullfile(wellviddir, [vidname '.mat']);
save(savename_matfile, 'LocoAVG','LocoSEM')
fprintf('-------------------------------------Finished-------------------------------------\n')

end