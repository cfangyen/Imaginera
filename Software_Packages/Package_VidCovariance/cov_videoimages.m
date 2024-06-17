function cov_vid = cov_videoimages(vid, sscale, tscale)
%COV_VIDEOIMAGES: Compute the covariance matrix of a sequence of images
sampleimg = im2gray(read(vid, 1));
sampleimg = imresize(sampleimg, sscale);
numFrames = vid.NumFrames;
height = size(sampleimg,1);
width = size(sampleimg,2);
T = 1:(1/tscale):numFrames;
allFrames_2d = zeros(numel(T), height*width);
for i=1:numel(T)
    currimg = im2gray(read(vid, T(i)));
    currimg = imresize(currimg, sscale);
    currimg = double(currimg);
    allFrames_2d(i,:) = currimg(:);
end
cov_vid = cosineSimilarity(allFrames_2d);
end

