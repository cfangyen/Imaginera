ParentDir = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol';
VidDir = dir(fullfile(ParentDir, '*.avi'));
NumVid = numel(VidDir);

% % Noise manipulation
k_noise = [0.0000, 0.0001, 0.0003, 0.0009, 0.0027, 0.0081, 0.0243, 0.0729];

% % Destination for noisy videos
NoisyDir = fullfile(ParentDir, 'Noisy videos');

for i = 1 : NumVid
    currVidName = VidDir(i).name;
    currNewName = sprintf('7NoisyVid%03d.avi', i); %%%
    currVid = VideoReader(fullfile(ParentDir, currVidName));
    FrameRate = currVid.FrameRate;
    AllFrames = read(currVid);
    if contains(currVid.VideoFormat, 'RGB')
        AllFrames = uint8(squeeze(mean(AllFrames, 3)));
    elseif contains(currVid.VideoFormat, 'Grayscale')
        AllFrames = unit8(AllFrames);
    end

    currOutVid = VideoWriter(fullfile(NoisyDir, currNewName), 'Uncompressed AVI');
    currOutVid.FrameRate = FrameRate;
    open(currOutVid);

    % Loop through each frame, add Gaussian, write to output video
    % pixmax = max(double(AllFrames(:)));
    % pixavg = mean(double(AllFrames(:)));
    % pixstd = std(double(AllFrames(:)));
    % srn = (pixmax - pixavg - pixstd)/(pixavg+pixstd);
    for j = 1 : size(AllFrames, 3)
        frame = AllFrames(:,:,j);
        
        % Add Gaussian noise with specified noise level
        noisyFrame = imnoise(frame, 'gaussian', 0, k_noise(end)); %%%
        
        % Write the noisy frame to the output video
        writeVideo(currOutVid, noisyFrame);
    end
    % Close the video writer object
    close(currOutVid);
end