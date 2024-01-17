p = 'C:\Users\fffei\Dropbox\Paper\Egg-laying\Data\2023-3-18 2hr 2min 15sec 20fps 6 well dop3\original';
savepath = 'C:\Users\fffei\Dropbox\Paper\Egg-laying\Data\WormTrainingData\DatasetImages';
viddir = dir(fullfile(p, '*.avi'));
load(fullfile(p, '..', 'masks6.mat'), 'roi')
numwells = numel(roi);
numvids = numel(viddir);
numdiv = 1;
for i = 1 : numvids
    vidname = fullfile(viddir(i).folder, viddir(i).name);
    [~, fname] = fileparts(vidname);
    vid = VideoReader(vidname);
    vidCol = vid.Width;
    vidRow = vid.Height;
    numFrames = vid.NumFrames;
    for k = 1 : 30 : numFrames
        currentImage = rgb2gray(read(vid, k));
        for j = 1 : numwells
            ROI = roi{j};
            angles = linspace(0, 2*pi, 10000);
            x = cos(angles) * ROI.Radius + ROI.Center(1);
            y = sin(angles) * ROI.Radius + ROI.Center(2);
            mask = poly2mask(x, y, vidRow, vidCol);
            maskedImage = currentImage;
            maskedImage(~mask) = nan;
            props = regionprops(mask, 'BoundingBox');
            maskedImage = imcrop(maskedImage, props.BoundingBox);
            [height, width] = size(maskedImage);
            sub_h = ceil(height / numdiv);
            sub_w = ceil(width / numdiv);
%             sub_images = cell(3,3);
            for row = 1 : numdiv
                for col = 1 : numdiv
                    r_start = (row - 1) * sub_h + 1;
                    r_end = mean(min([height r_start+sub_h-1]));
                    c_start = (col - 1) * sub_w + 1;
                    c_end = mean(min([width c_start+sub_w-1]));
                    
                    sub_image = maskedImage(r_start:r_end, c_start:c_end);
%                     sub_images{row, col} = sub_image;
                    date = datetime;
                    imwrite(sub_image, fullfile(savepath,...
                        [sprintf('3-20_R%d_C%d_W%d_F%d', row,col,j,k) fname '.png']))
                end
            end
        end
    end
end