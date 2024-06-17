function vid2d_pca_removed = PCA_bkgsubtraction(vidpath, resize_ratio, do_savevid)
% Read video
[path2save,oldvidname,ext] = fileparts(vidpath);
vid = VideoReader(vidpath);
sampleimg = read(vid, 1);
sampleimg_resize = imresize(sampleimg, resize_ratio);
numFrames = vid.NumFrames;
height = size(sampleimg_resize, 1);
width = size(sampleimg_resize, 2);
duration = vid.Duration;
fps = numFrames/duration;
allFrames_2d = zeros(height*width,numFrames);
fprintf('Start reading video: %s\n', oldvidname);
for i=1:numFrames
    if mod(i, 1000) == 0
        fprintf('Read %d out of %d\n', i, numFrames)
    end
    currimg = imresize(read(vid, i), resize_ratio);
    allFrames_2d(:,i) = currimg(:);
end
fprintf('Finish reading video\n');
% Perform PCA on the 2D matrix
fprintf('Start PCA of video: %s\n', oldvidname);
[coeff,score,~] = pca(allFrames_2d','Centered',false);
fprintf('PCA finished\n');
% The first principal component
first_pc = coeff(:,1);
% Subtract the first principal component from the original 2D matrix
allFrames_2d = allFrames_2d - first_pc * score(:,1)';
vid2d_pca_removed = allFrames_2d;
% Reshape the 2D matrix back to the original 3D array size
allFrames_3d = reshape(allFrames_2d, [height, width, numFrames]);
if do_savevid
    % write video
    newvidname = [oldvidname '_pca' ext];
    fprintf('Saving video: %s\n', newvidname);
    allFrames_3d = uint8(allFrames_3d);
    vidpath2save = fullfile(path2save, newvidname);
    v = VideoWriter(vidpath2save);
    v.FrameRate = fps;
    open(v);
    for i = 1 : numFrames
        frame = im2gray(allFrames_3d(:,:,i));
        writeVideo(v,frame);
    end
    close(v)
    fprintf('Video saved\n');
end
end

