function I = Hu_calculation_mw(WellVid_3d, SScale, TScale)
%% Initalizing parameters
startIndex = 1;
endIndex = size(WellVid_3d, 3);
numFrames = endIndex - startIndex + 1;
StepFrames = round(1/TScale);
sampleFrames = 1 : StepFrames : numFrames;
sigma = 1;
numObs = numel(sampleFrames);
I = zeros(numObs, 8);

okay = 1;
while okay == 0
    Img = WellVid_3d(:,:, startIndex);
    smoothedImg = imgaussfilt(Img, sigma);
    figure;
    imagesc(smoothedImg); colormap("gray")
    set(gcf, 'Position', [232,246,500,500])

    okay = input('If sigma and threshold okay, enter 1, otherwise enter 0:  ');

    if okay==0
        sigma = input('Enter new sigma: ');
    end
    close all
end

% define a coordinate system for image
% sampleimg = im2gray(read(WellVid_3d,1));
sampleimg = WellVid_3d(:,:, 1);
sampleimg = imresize(sampleimg, SScale, "bilinear");
[height, width] = size(sampleimg);
xgrid = repmat((-floor(height/2):1:ceil(height/2)-1)',1,width);
ygrid = repmat(-floor(width/2):1:ceil(width/2)-1,height,1);

for i = 1 : numObs
    j = (i-1)*StepFrames + startIndex;

    currentimg = double(WellVid_3d(:,:, j));
    % Reduce noise
    smoothedimg = currentimg;
    if SScale < 1
        currentimg_sm = imresize(smoothedimg, SScale, "bilinear");
    else
        currentimg_sm = smoothedimg;
    end
    cm_00 = sum(currentimg_sm, 'all');

    [x_bar, y_bar] = centerOfMass(currentimg_sm,cm_00,xgrid,ygrid);
    % normalize coordinate system by subtracting mean
    xnorm = x_bar - xgrid;
    ynorm = y_bar - ygrid;
    % central moments
    xnorm0 = xnorm.^0;      ynorm0 = ynorm.^0;
    xnorm1 = xnorm0.*xnorm; ynorm1 = ynorm0.*ynorm;
    xnorm2 = xnorm1.*xnorm; ynorm2 = ynorm1.*ynorm;
    xnorm3 = xnorm2.*xnorm; ynorm3 = ynorm2.*ynorm;

    mu_11 = central_moments(currentimg_sm,cm_00,xnorm1,ynorm1,1,1);
    mu_20 = central_moments(currentimg_sm,cm_00,xnorm2,ynorm0,2,0);
    mu_02 = central_moments(currentimg_sm,cm_00,xnorm0,ynorm2,0,2);
    mu_21 = central_moments(currentimg_sm,cm_00,xnorm2,ynorm1,2,1);
    mu_12 = central_moments(currentimg_sm,cm_00,xnorm1,ynorm2,1,2);
    mu_03 = central_moments(currentimg_sm,cm_00,xnorm0,ynorm3,0,3);
    mu_30 = central_moments(currentimg_sm,cm_00,xnorm3,ynorm0,3,0);
    % Hu's Invariant moments
    I(i,1) = mu_20 + mu_02;
    I(i,2) = (mu_20  - mu_02 )^2 + 4*(mu_11 )^2;
    I(i,3) = (mu_30  - 3*mu_12 )^2 + (mu_03  - 3*mu_21 )^2;
    I(i,4) = (mu_30  + mu_12 )^2 + (mu_03  + mu_21 )^2;
    I(i,5) = (mu_30  - 3*mu_12 )*(mu_30  + mu_12 )*((mu_30  + mu_12 )^2 - 3*(mu_21  + mu_03 )^2) + (3*mu_21  - mu_03 )*(mu_21  + mu_03 )*(3*(mu_30  + mu_12 )^2 - (mu_03  + mu_21 )^2);
    I(i,6) = (mu_20  - mu_02 )*((mu_30  + mu_12 )^2 - (mu_21  + mu_03 )^2) + 4*mu_11 *(mu_30  + mu_12 )*(mu_21  + mu_03 );
    I(i,7) = (3*mu_21  - mu_03 )*(mu_30  + mu_12 )*((mu_30  + mu_12 )^2 - 3*(mu_21  + mu_03 )^2) + (mu_30  - 3*mu_12 )*(mu_21  + mu_03 )*(3*(mu_30  + mu_12 )^2 - (mu_03  + mu_21 )^2);
    I(i,8) = mu_11 *(mu_30  + mu_12 )^2 - (mu_03  + mu_21 )^2 - (mu_20  - mu_02 )*(mu_30  + mu_12 )*(mu_21  + mu_03 );
end

end


function [x_bar, y_bar] = centerOfMass(image,cm_00,xgrid,ygrid)

eps = 10^(-6); % very small constant

x_bar = sum((xgrid.*image), 'all')/(cm_00+eps);
y_bar = sum((ygrid.*image), 'all')/(cm_00+eps);

end

function cm = central_moments(image,cm_00,xnormp,ynormq,p,q)

cm = sum((xnormp).*(ynormq).*image, 'all');
% cm_00 = sum(image, 'all'); %this is same as mu(0,0);
% normalise moments for scale invariance
cm = cm/(cm_00^(1+(p+q)/2));

end
