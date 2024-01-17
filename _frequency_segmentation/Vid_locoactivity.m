function [LocoAVG, LocoSEM] = Vid_locoactivity(Vidname, separation, background, roi)
vidObj = VideoReader(Vidname);
FrameSample = rgb2gray(read(vidObj,1));
[numRow, numCol] = size(FrameSample);
[~, Vname, ~] = fileparts(Vidname);
numFrames = vidObj.NumFrames;
FrameRate = vidObj.FrameRate;
numwells = numel(roi);
LocoAVG = nan(numwells, numel(separation));
LocoSEM = nan(numwells, numel(separation));

angles = linspace(0, 2*pi, 10000);

% Read all frames from the video and store them in memory
allFrames = zeros(numRow, numCol, numFrames, 'uint8');
for frameIdx = 1 : numFrames
    allFrames(:,:,frameIdx) = rgb2gray(read(vidObj, frameIdx)) - background;
end
parfor j = 1 : numel(separation)
    currentSep = separation(j);
    LocoAVG_cum_currentsep = zeros(numwells, numFrames-currentSep);
    for i = 1 : numFrames-currentSep
        previousImage = allFrames(:,:,i);
        currentImage  = allFrames(:,:,min([numFrames i+currentSep]));

        for kk = 1 : numwells
            fprintf('Well = %d, SEP = %d, Analyzing video %s, frame%d out of %d\n',kk, currentSep, Vname, i, numFrames-currentSep)
            ROI = roi{kk};
            x = cos(angles) * ROI.Radius + ROI.Center(1);
            y = sin(angles) * ROI.Radius + ROI.Center(2);
            mask = poly2mask(x, y, numRow, numCol);
            props = regionprops(mask, 'BoundingBox');
            
            % crop images with mask
            maskedPrev = previousImage;
            maskedPrev(~mask) = nan;
            maskedPrev = imcrop(maskedPrev, props.BoundingBox);
            
            maskedcurr = currentImage;
            maskedcurr(~mask) = nan;
            maskedcurr = imcrop(maskedcurr, props.BoundingBox);
            
            pixeldiff = abs(maskedPrev - maskedcurr);
            pxlmeanint = mean([maskedPrev(:);maskedcurr(:)], 'omitnan');
            pixeldiff_norm = pixeldiff/pxlmeanint;
            diffImg = pixeldiff_norm;
            diffImg_gauss = imgaussfilt(diffImg, 1);
            T = graythresh(diffImg_gauss);
            diffBW = imbinarize(diffImg_gauss, T);
            diffImg_gauss(~diffBW) = nan;
            DIFF = sum(diffImg_gauss(:), 'omitnan');
            LocoAVG_cum_currentsep(kk, i) = DIFF;
        end
    end
    LocoAVG_cum_currentsep_persec = LocoAVG_cum_currentsep*FrameRate;
    LocoAVG(:, j) = mean(LocoAVG_cum_currentsep_persec, 2);
    LocoSEM(:, j) = std(LocoAVG_cum_currentsep_persec,0,2)/sqrt(numFrames-currentSep);
end
delete(gcp('nocreate'))
end

