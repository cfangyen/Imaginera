clear;clc;
close all;
%% video initialization
vidpath='D:\dropbox\Dropbox\Dian Hongfei shared\worm analysis\Data\High-Res Swim\Exp_MT14984\MT149840001_bkgsubtracted.avi';
vid= VideoReader(vidpath);
numFrames=vid.NumFrames;
height=vid.Height;
width=vid.Width;
duration=vid.Duration;
fps=round(numFrames/duration);
allFrames_3d=zeros(height,width,numFrames);
for i=1:numFrames
    allFrames_3d(:,:,i)=vid.readFrame;
end
%% PCA
allFrames_2d = reshape(allFrames_3d, [height*width, numFrames]);

% Perform PCA on the 2D matrix
[coeff,score,latent]=pca(allFrames_2d','Centered',false);

% The first principal component
first_pc = coeff(:,1);

% Reshape the first principal component back to the original image size
first_pc_image = reshape(first_pc, [height, width]);

% Display the first principal component
figure; imshow(first_pc_image, []);
% Subtract the first principal component from the original 2D matrix
allFrames_2d = allFrames_2d - first_pc * score(:,1)';
% Reshape the 2D matrix back to the original 3D array size
allFrames_3d = reshape(allFrames_2d, [height, width, numFrames]);
%% write video
allFrames_3d=uint8(allFrames_3d);
vidpath='D:\dropbox\Dropbox\Dian Hongfei shared\worm analysis\Data\High-Res Swim\Exp_MT14984\MT149840001_bkgsubtracted_PCA1.avi';
v=VideoWriter(vidpath);
v.FrameRate = fps;
open(v);

for i=1:numFrames
    frame=im2gray(allFrames_3d(:,:,i));
    writeVideo(v,frame);
end
close(v)