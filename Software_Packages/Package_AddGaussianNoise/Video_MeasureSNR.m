ParentDir = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol\Noisy videos';
% % Noise manipulation
k_noise = [0.0001, 0.0003, 0.0009, 0.0027, 0.0081, 0.0243];
noiseDir = {'p0000'};
snrValues = zeros(size(k_noise));
wormsize = 10;
for j = 1 : numel(k_noise)
    currDir = fullfile(ParentDir, noiseDir{1});
    VidDir = dir(fullfile(currDir, '*.avi'));
    NumVid = numel(VidDir);
    curr_snr = zeros(1, 1);
    for i = 1 %: NumVid
        currVidName = VidDir(i).name;
        fprintf('Noise level: %.4f, Video %03d\n', k_noise(j), i)
        currVid = VideoReader(fullfile(currDir, currVidName));
        AllFrames = read(currVid);
        if contains(currVid.VideoFormat, 'RGB')
            AllFrames = uint8(squeeze(mean(AllFrames,3)));
        elseif contains(currVid.VideoFormat, 'Grayscale')
            AllFrames = uint8(AllFrames);
        end
        signalPower  = mean(((max(double(AllFrames),[],[1 2]))).^2)/wormsize;
        AllFrames_noisy = zeros(size(AllFrames));
        for k = 1 : size(AllFrames, 3)
            frame = AllFrames(:,:,k);
            noisyFrame = imnoise(frame, 'gaussian', 0, k_noise(j)); %%%
            AllFrames_noisy(:,:,k) = noisyFrame;
        end
        noisePower = var(double(AllFrames_noisy(:))-double(AllFrames(:)));
        curr_snr(i) = 10*log10(signalPower / noisePower);
    end
    snrValues(j) = mean(curr_snr);
end
figure(1);clf
plot(k_noise, snrValues)
set(gca, 'XScale', 'log')
%% Adding noise to clean image
% Load the clean image
cleanImagePath = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\Manuscript\Individual Figures\Noise reduced 0.03 s.tif';
cleanImage = imread(cleanImagePath);
% Load background noise
bkgPath = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\Manuscript\Individual Figures\noise.bmp';
noiseImage = imread(bkgPath);

% Add noise to the clean image
noisyImage = double(cleanImage) + double(noiseImage);

% Clip the noisy image to valid range
noisyImage = uint8(min(max(noisyImage, 0), 255));

% Display the images
figure(1); clf
imshow(noisyImage, []);
title('Noisy Image');

% Save the noisy image
imwrite(noisyImage, 'noisy_version_of_image.tif');
disp('Noisy image created and saved successfully.');

