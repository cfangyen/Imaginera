function J = MSER_calculation_mw(WellVid_3d, SScale, TScale,ThresholdDelta,MaxAreaVariation, RegionAreaRange, do_plot)
%% Initalizing parameters
startIndex = 1;
endIndex = size(WellVid_3d, 3);
numFrames = endIndex - startIndex + 1;
StepFrames = round(1/TScale);
sampleFrames = 1 : StepFrames : numFrames;
numObs = numel(sampleFrames);
J = cell(numObs, 1);

for i = 1 : numObs
    j = (i-1)*StepFrames + startIndex;
%     currentimg = im2gray(read(WellVid_3d,j));
    currentimg = WellVid_3d(:,:, j);
    if SScale < 1
        currentimg = imresize(currentimg, SScale, "bilinear");
    end
%     points = detectORBFeatures(currentimg,'ScaleFactor',1.01,'NumLevels',50);
    [regions, ~] = detectMSERFeatures(currentimg, "ThresholdDelta", ThresholdDelta, 'MaxAreaVariation', MaxAreaVariation, "RegionAreaRange", RegionAreaRange);
    % test plot
    if do_plot
        %%
        figure(1); clf
        imagesc(currentimg);
        hold on;
        plot(regions, 'showPixelList',true,'showEllipses',true);
        %     plot(points);
        hold off;
        pause(0.001)
        %% to delete
        figure(2); clf
        imagesc(currentimg); axis image; colormap('turbo')
        hold on
        plot(regions, 'showPixelList',false,'showEllipses',true);
        xticks([])
        yticks([])
        xlim([120 175])
        ylim([210 280])
%         xlabel('X position')
%         ylabel('Y position')
        set(gca, 'FontName', 'Arial', 'FontSize', 6)
        set(gcf, 'Position', [180,163,400,600])
        %% to delete
    end
    J{i} = regions;
    % % Write this frame out to a new video file.
    % writeVideo(outputVideo, getframe(gcf));
end
% % Close video writer
% close(outputVideo);
end

