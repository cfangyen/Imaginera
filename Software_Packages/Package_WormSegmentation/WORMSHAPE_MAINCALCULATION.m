function curvdatafiltered = WORMSHAPE_MAINCALCULATION(vidObj, options, resizefactor, tscale)
%WORMSHAPE_MAINCALCULATION Analyze the undulatory dynamics of worms
wormlabel     = options{2};
fps           = options{3}; fps = fps*tscale;
pix_per_mm    = options{4}; %
wormthreshold = options{5};
thisperiod    = options{6};
decim         = options{7};
filsize       = options{8}; %
spline_p      = options{11};
expo = 0.7;
istart = thisperiod(1);
iend   = thisperiod(2);
skip = floor((iend-istart+1)/10);

invert_img   = 0;
decim_filter = ones(decim) / (decim^2);

if skip ==0
    skip = 1;
end
idxskip = round(1/tscale);
% numframes  = iend - istart + 1;
numframes = numel(istart: idxskip: iend);

numcurvpts = 100;

mov_size_multiplier = 1;
savefps = 30;
mov_quality = .9;



%%%%%%%% PREVIEW IMAGES %%%%%%%%%%%
%%% Worm Analysis
j=0;
% manually remove bright none-worm objects
img = mean(read(vidObj,istart),3);
img = imresize(img, resizefactor, 'bilinear');
img = imfilter(img, decim_filter, 'same');
img = img(1:decim:end,1:decim:end);
bw_remove = false(size(img));
%%%%%%%%%%%%%%%%%%%%%%%
img(bw_remove) =  min(min(img));
% % manually select four pixel points from the background to eliminate the
% % background noise
% figure(2); clf;
% imagesc(img); colormap gray; axis image; hold on;
% title('Pick four points as background pixels')
% [xs_bkg, ys_bkg] = ginput(4);
% plot(xs_bkg, ys_bkg, 'or')

% create sum image
for i=istart:skip:iend
    j = j+1;
    
    img = mean(read(vidObj,i),3);
    img = imresize(img, resizefactor, 'bilinear');
    
    img = imfilter(img, decim_filter, 'same');
    img = img(1:decim:end,1:decim:end);
    
    if invert_img
        img = 255-img;
    end
    if i == istart
        imgsum = single(img);
        [ysize, xsize] = size(img);
        imgmin = ones(size(img));
        imgdata = zeros(ysize, xsize, length(istart:skip:iend));
    end
    % figure(1);
    % imagesc(img); colormap gray;hold on;
    % axis image;    title(num2str(i));
    imgdata(:,:,j) = img;
    imgsum = imgsum + single(img);
end


% figure(1);clf;
% imagesc(imgsum); colormap jet; hold on;
% title('sum image');

text(10,20, 'select ROI: upper left then lower right', 'Color', 'white');
% [cropx1, cropy1] = ginput(1);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to release
    cropx1 = 1; cropy1 = 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cropx1 = floor(cropx1);
cropy1 = floor(cropy1);

% ADF EDIT: Make sure the crops are in-bounds.
cropx1 = max([1,cropx1]);
cropy1 = max([1,cropy1]);
cropx1 = min([size(img,2),cropx1]);
cropy1 = min([size(img,1),cropy1]);

% Get the second corner of the ROI
plot([1 xsize], [cropy1 cropy1], '-r');
plot([cropx1 cropx1], [1 ysize], '-r');
% [cropx2, cropy2 ] = ginput(1);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to release
cropx2 = size(img,2); cropy2 = size(img,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cropx2 = floor(cropx2);
cropy2 = floor(cropy2);

% ADF EDIT: Make sure the crops are in-bounds
cropx2 = max([1,cropx2]);
cropy2 = max([1,cropy2]);
cropx2 = min([size(img,2),cropx2]);
cropy2 = min([size(img,1),cropy2]);

plot([1 xsize], [cropy2 cropy2], '-r');
plot([cropx2 cropx2], [1 ysize], '-r');

%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter Initiation %%%%%%%%%%%%%%%%%%%%%%%%%%%

showcalc = 0;
deinterlace = 1; interlaceframe = 1;
cropyes = 1;


curvdata = nan(numframes,numcurvpts);
areadata = zeros(numframes,1);
centroiddata =  zeros(numframes,2);

cv2i_data = zeros(numframes,numcurvpts+2,2);
angledata = zeros(numframes,numcurvpts+1);
path1_rescaled_data = zeros(numframes,numcurvpts,2);
path2_rescaled_data = zeros(numframes,numcurvpts,2);
lendata = zeros(numframes, 1);

for j=1
    i = istart + (j - 1);
    img = mean(read(vidObj,i),3);%%%%%%%%%%%%%%%% changed
    img = imresize(img, resizefactor, 'bilinear');
    
    img = imfilter(img, decim_filter, 'same');
    img = img(1:decim:end,1:decim:end);
    if invert_img
        img = 255-img;
    end
    
    img2 = abs(single(img(:,:,1))- imgmin);
    img = abs(single(img(:,:,1)));
    
    if deinterlace
        img(3-interlaceframe:2:end) = img(interlaceframe:2:end);
        img2(3-interlaceframe:2:end) = img2(interlaceframe:2:end);
    end
    if cropyes
        imgcrop = img(cropy1:cropy2,cropx1:cropx2);%
        
        imgcrop2 = img2(cropy1:cropy2,cropx1:cropx2);%
        
    else
        imgcrop = img;
        imgcrop2 = img2;
    end
    imgcrop = imgcrop';
    imgcrop2 = imgcrop2';
    
    imgcrop3 = imgcrop - imgcrop2;
    imgcrop4 = imgcrop3 - min(min(imgcrop3));
    [a,c] = find (imgcrop4 > 40);
    contour = [a,c];
    % figure(2); hold off;
    % imagesc(imgcrop4, [5 250]); hold on;
    % colormap gray
end

ddd1 = [];
vvv1 = [];

for j=1:numframes
    i = istart + (j - 1)*idxskip;
    
    if i>vidObj.NumberOfFrames; break; end
%     
    img = mean(read(vidObj,i),3);
    img(bw_remove) =  min(min(img));
    img = imresize(img, resizefactor, 'bilinear');
    img = imfilter(img, decim_filter, 'same');
    img = img(1:decim:end,1:decim:end);
    if invert_img
        img = 255-img;
    end
    
    img = abs(single(img(:,:,1))- imgmin);
    img2 = abs(single(img(:,:,1)));
    
    if deinterlace
        img(3-interlaceframe:2:end) = img(interlaceframe:2:end);
        img2(3-interlaceframe:2:end) = img2(interlaceframe:2:end);
        
    end
    if cropyes
        imgcrop = img(cropy1:cropy2,cropx1:cropx2);
        imgcrop2 = img2(cropy1:cropy2,cropx1:cropx2);
    else
        imgcrop = img;
        imgcrop2 = img2;
    end
    imgcrop = imgcrop';
    imgcrop2 = imgcrop2';
    
    
    if j==1
        % colormap jet;
        [ysize, xsize ] = size(imgcrop);
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to remove
        headx = 1; heady = 1;
        headx0 = headx; heady0 = heady;
        tailx = xsize; taily = ysize;
        tailx0 = tailx; taily0 = taily;

        
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to release
        worm_diam = 5*resizefactor;
        % title(['worm diameter = ' num2str(worm_diam) ' pixels']);
        worm_area_est = 10*worm_diam^2;
        sizethresh = round(worm_area_est / 2);
        if mod(round(filsize*worm_diam),2)==1
            filradius = round(filsize*worm_diam/2);
        else
            filradius = round(filsize*worm_diam/2)+1;
        end

        if filradius==0
            filradius = 1;
        end
        
        fil = fspecial('disk', filradius);

        colormap gray;
        zoom out;
        set(gcf, 'Position',    [ 129   190   310   463]);
        
        
    end
    
    img2 = conv2(single(imgcrop), fil, 'same');
    
    lvl = min(min(img2))+wormthreshold* (-min(min(img2))+max(max(img2)));
    
    bw =(img2> lvl);
    
    
    bw2 = bwareaopen(bw,  sizethresh);
    bw3 = imcomplement(bw2);
    bw4 = bwareaopen(bw3, sizethresh);
    bw5 = imcomplement(bw4);
    
    STATS = regionprops(logical(bw5),'Area', 'Centroid');
    
    if size(STATS,1) == 0
        disp('Error: no worm found');
        continue;
    end
    
    areadata(j) = STATS.Area;
    centroiddata(j,:) = STATS.Centroid;
    B = bwboundaries(bw5, 'noholes'); %  trace boundary clockwise
    
    B1 = B{1}; % boundary coordinates
    
    B1_size = size(B1,1);
    
    ksep = ceil(B1_size/20);
    
    B1_plus = circshift(B1,[ksep 0]);
    B1_minus = circshift(B1,[-ksep 0]);
    
    AA = B1 - B1_plus;  % AA and BB are vectors between a point on boundary and neighbors +- ksep away
    BB = B1 - B1_minus;
    
    cAA = AA(:,1) + sqrt(-1)*AA(:,2);
    cBB = BB(:,1) + sqrt(-1)*BB(:,2);
    
    B1_angle = unwrap(angle(cBB ./ cAA));
    
    min1 = find(B1_angle == min(B1_angle),1); % find point on boundary w/ minimum angle between AA, BB
    B1_angle2 = circshift(B1_angle, -min1);
    min2a = round(.25*B1_size)-1+find(B1_angle2(round(.25*B1_size):round(0.75*B1_size))==min(B1_angle2(round(.25*B1_size):round(0.75*B1_size))),1);  % find minimum in other half
    min2 = 1+mod(min2a + min1-1, B1_size);
    
    tmp = circshift(B1, [1-min1 0]);
    end1 = 1+mod(min2 - min1-1, B1_size);
    
    path1 = tmp(1:end1,:);
    path2 = tmp(end:-1:end1,:);
    
    if norm(path1(1,:) - [heady headx]) > norm(path1(end,:) - [heady headx]) % if min1 is at tail, reverse both paths
        tmp = path1;
        path1 = path2(end:-1:1,:);
        path2 = tmp(end:-1:1,:);
    end
    
    heady = path1(1,1);
    headx = path1(1,2);
    taily = path1(end,1);
    tailx = path1(end,2);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to release
    
    path_length = numcurvpts;
    
    path1_rescaled = zeros(path_length,2);
    path2_rescaled = zeros(path_length,2);
    path1_rescaled2 = zeros(path_length,2);
    path2_rescaled2 = zeros(path_length,2);
    
    path1_rescaled(:,1) = interp1(0:size(path1,1)-1, path1(:,1), (size(path1,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path1_rescaled(:,2) = interp1(0:size(path1,1)-1, path1(:,2), (size(path1,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path2_rescaled(:,1) = interp1(0:size(path2,1)-1, path2(:,1), (size(path2,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path2_rescaled(:,2) = interp1(0:size(path2,1)-1, path2(:,2), (size(path2,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    
    
    for kk=1:path_length
        tmp1 = repmat(path1_rescaled(kk,:), [path_length,1]) - path2_rescaled;
        tmp2 = sqrt(tmp1(:,1).^2 + tmp1(:,2).^2);
        path2_rescaled2(kk,:) = path2_rescaled(find(tmp2==min(tmp2),1),:);
    end
    
    for kk=1:path_length
        tmp1 = repmat(path2_rescaled(kk,:), [path_length,1]) - path1_rescaled;
        tmp2 = sqrt(tmp1(:,1).^2 + tmp1(:,2).^2);
        path1_rescaled2(kk,:) = path1_rescaled(find(tmp2==min(tmp2),1),:);
    end
    
    
    
    dorsalx = path1_rescaled2(:,1);
    dorsaly = path1_rescaled2(:,2);
    ventralx = path2_rescaled2(:,1);
    ventraly = path2_rescaled2(:,2);
    
    dorsal = [ventralx,ventraly];
    ventral = [dorsalx,dorsaly];
    
    dorsalline =  round(dorsal);
    ventralline = round(ventral);
    
    
    a2 =[];
    a3 =[];
    
    
    for i = 1:length(ventralline)
        a1 = find(ventralline(i,1) == contour(:,1) & ventralline(i,2) == contour(:,2));
        a4 = find(dorsalline(i,1) == contour(:,1) & dorsalline(i,2) == contour(:,2));
        if a1 > 0
            a1 = 1;
        else
            a1 = 0;
        end
        a2 = cat(2,a2,a1);
        
        if a4 > 0
            a4 = 1;
        else
            a4 = 0;
        end
        a3 = cat(2,a3,a4);
        
    end
    ddd = sum(a3);
    vvv = sum(a2);
    ddd1 = cat(1,ddd1,ddd);
    vvv1 = cat(1,vvv1,vvv);
    comb_cont = cat(1,ddd1',vvv1');
    
    weight_fn = ones(path_length,1);
    tmp=round(path_length*0.2);
    weight_fn(1:tmp)=(0:tmp-1)/tmp;
    weight_fn(end-tmp+1:end)=(tmp-1:-1:0)/tmp;
    weight_fn = [weight_fn weight_fn];
    
    midline   = 0.5*(path1_rescaled+path2_rescaled);
    midline2a = 0.5*(path1_rescaled+path2_rescaled2);
    midline2b = 0.5*(path1_rescaled2+path2_rescaled);
    midline_mixed = midline2a .* weight_fn + midline .* (1-weight_fn);
    
    Line = midline_mixed;
    
    interpfactor = 10;
    Line2 = interp1(Line, (1:(1/interpfactor):100)); % worm's center line in xy coordinates in imgcrop
    
    xy = circshift(Line2, [0 1])'; df = diff(xy,1,2);
    
    
    
    t = cumsum([0, sqrt([1 1]*(df.*df))]);
    cv = csaps(t,xy,spline_p);
    
    dorsal_xy = circshift(dorsalline, [0 1])'; df = diff(dorsal_xy,1,2);
    if ~any(df(:))
        continue
    end
    tmpt2 = cumsum([0, sqrt([1 1]*(df.*df))]);
    dorsal_cv = csaps(tmpt2,dorsal_xy,spline_p);
    
    ventral_xy = circshift(ventralline, [0 1])'; df = diff(ventral_xy,1,2);
    if ~any(df(:))
        continue
    end
    tmpt2 = cumsum([0, sqrt([1 1]*(df.*df))]);
    ventral_cv = csaps(tmpt2,ventral_xy,spline_p);

    if showcalc
        figure(2);
        imshow(bw); hold on
        plot(path1_rescaled(1,2),  path1_rescaled(1,1), 'or');
        plot(path2_rescaled(end,2),path2_rescaled(end,1), 'og');
        fnplt(cv, '-g');
        plot(xy(1,end),xy(2,end),'*r'); % centre point of cv
        plot(xy(1,1),xy(2,1),'*b');
        fnplt(dorsal_cv,  '-r');
        fnplt(ventral_cv, '-r');
        hold off
        axis tight
        set(gcf, 'Position', [235,273,554,447])
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Original image
%     figure(2); clf; axis image
%     hold on
%     J = imnoise(im2double(imgcrop)/255*1.8, 'gaussian',.2,0.005);
%     imshow(J)
%     colormap("gray")
%     xlim([190 330])
%     ylim([270 455])
%     set(gcf, 'Position', [345,649,200,250])
    %%%%%%%%%%%%%%%%%%%%%%%%% Image after PCA
%     figure(3); clf; axis image
% %     hold on
%     imshow(imgcrop'/255*2)
%     colormap("gray")
%     xlim([260 460]*resizefactor)
%     ylim([170 370]*resizefactor)
%     drawnow
%     set(gcf, 'Position', [345,649,200,250])
    % pause(0.1)
%     %%%%%%%%%%%%%%%%%%%%%%%%% Thresholded image 
%     figure(4); clf; axis image
%     imshow(img2> lvl*1.5);
%     hold on
%     fnplt(cv, '-g');
%     plot(xy(1,end),xy(2,end),'*r'); % centre point of cv
%     plot(xy(1,1),xy(2,1),'*b');
%     fnplt(dorsal_cv,  '--r');
%     fnplt(ventral_cv, '--r');
%     hold off
%     xlim([190 330])
%     ylim([270 455])
%     set(gcf, 'Position', [345,649,200,250])
%     %%%%%%%%%%%%%%%%%%%%%%%%% Worm segmentation
%     figure(5); clf; axis image
%     imshow(img2> 255);
%     hold on
%     fnplt(cv, '-g');
%     plot(xy(1,end),xy(2,end),'*r'); % centre point of cv
%     plot(xy(1,1),xy(2,1),'*b');
%     fnplt(dorsal_cv,  '--r');
%     fnplt(ventral_cv, '--r');
%     hold off
%     xlim([190 330])
%     ylim([270 455])
%     set(gcf, 'Position', [345,649,200,250])
%     %%%%%%%%%%%%%%%%%%%%%%%%% surf plot of worm image
%     figure(6); clf; axis image
%     imagesc(img2);
%     colormap("turbo")
%     xticks([200 250 300])
%     yticks([300 350 400 450])
%     xlabel('X position')
%     ylabel('Y position')
%     xlim([190 330])
%     ylim([270 455])
%     set(gcf, 'Position', [345,649,144,180])
%     set(gca, 'FontName', 'Arial', 'FontSize', 6)
    %%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     if domovie && j>1
%         MakeQTMovie('addframe');
%     end
%     
%     if j==1
%         plot([Line(1,2) headx0],[Line(1,1) heady0], '-oc');
%         plot([Line(end,2) tailx0],[Line(end,1) taily0], '-oc');
%         pause(1);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to release
    
    cv2 =  fnval(cv, t)';
    df2 = diff(cv2,1,1); df2p = df2';
    
    splen = cumsum([0, sqrt([1 1]*(df2p.*df2p))]);
    lendata(j) = splen(end);
    % interpolate to equally spaced length units
    cv2i = interp1(splen+.00001*(0:length(splen)-1),cv2, (0:(splen(end)-1)/(numcurvpts+1):(splen(end)-1)));

    % store cv2i data
    
    cv2i_data(j,:,:) = cv2i;
    path1_rescaled_data(j,:,:) = path1_rescaled;
    path2_rescaled_data(j,:,:) = path2_rescaled;
    df2 = diff(cv2i,1,1);
    atdf2 =  unwrap(atan2(-df2(:,2), df2(:,1)));
    curv = unwrap(diff(atdf2,1));
    curvdata(j,:) = curv';
    % calculate the angle of attack during worm's locomotion
    atdf2 = atan2(-df2(:,2), df2(:,1));
    theta   = (mean(max(atdf2)) + mean(min(atdf2)))/2;
    xcenter = cv2i(1,1);
    ycenter = cv2i(1,2);
    center  = repmat([xcenter ycenter], size(cv2i, 1), 1);
    Ro      = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    % do the rotation
    cv2io = (Ro*(cv2i' - center') + center')';
    df2o   = diff(cv2io,1,1);
    atdf2o = atan2(-df2o(:,2), df2o(:,1));
    angledata(j,:) = atdf2o';
end  % end main loop

% Post-processing raw curvature data
lenavg = mean(lendata);
curvdata = curvdata.*lenavg;
timefilter = max(1, round(fps*1/6));
bodyfilter = 5;
curvdata_median = medfilt2(curvdata, [timefilter bodyfilter]);
tmp = reshape(curvdata_median, [numel(curvdata_median),1]);
curv05 = prctile(tmp, 5);
curv95 = prctile(tmp, 95);
curvdata(curvdata > curv95) = curvdata_median(curvdata > curv95);
curvdata(curvdata < curv05) = curvdata_median(curvdata < curv05);
curvfilter = fspecial('average',[timefilter,bodyfilter]);
curvdatafiltered = imfilter(curvdata, curvfilter, 'replicate');
dKdt_data        = gradient(curvdatafiltered')'*fps;
% Post-processing raw angle data
angledata_median = medfilt2(angledata, [timefilter bodyfilter]);
tmp = reshape(angledata_median, [numel(angledata_median),1]);
angle05 = prctile(tmp, 5);
angle95 = prctile(tmp, 95);
angledata(angledata > angle95) = angledata_median(angledata > angle95);
angledata(angledata < angle05) = angledata_median(angledata < angle05);
anglefilter = fspecial('average',[timefilter,bodyfilter]);
angledatafiltered = imfilter(angledata,  anglefilter , 'replicate');
w_diam = worm_diam/pix_per_mm;


end

