function J = MSER_calculation(Vid3d, SScale, TScale, ThresholdDelta, MaxAreaVariation, RegionAreaRange, do_plot)
%% Initalizing parameters
startIndex = 1;
endIndex = size(Vid3d, 3);
numFrames = endIndex - startIndex + 1;
StepFrames = round(1/TScale);
sampleFrames = 1 : StepFrames : numFrames;
numObs = numel(sampleFrames);
J = cell(numObs, 1);

for i = 1 : numObs
    j = (i-1)*StepFrames + startIndex;
    currentimg = Vid3d(:,:,j);
    if SScale<1
        currentimg = imresize(currentimg, SScale, 'bilinear');
    end
%     points = detectORBFeatures(currentimg,'ScaleFactor',1.01,'NumLevels',50);
    [regions, ~] = detectMSERFeatures(currentimg', "ThresholdDelta", ThresholdDelta, 'MaxAreaVariation', MaxAreaVariation, "RegionAreaRange", RegionAreaRange);
    % test plot
    if do_plot %&& (i==round(17.7*30)||i==round(17.87*30)||i==round(18.1*30)||i==(18.3*30)||i==round(18.5*30)||i==round(18.7*30)||i==round(18.87*30))
        %%
        figure(1); 
        clf
        imagesc(currentimg');axis image;colormap("turbo")
        hold on;
        if ~isempty(regions)
            plot(regions(floor(linspace(1,regions.Count,3))), 'showPixelList',false,'showEllipses',true);
        end
        hold off;
        
%         xlim([20 90])
        % ylim([40 160])
        xticks([])
        yticks([])
        set(gcf, 'Position', [345,349,500,500])
        set(gca, 'FontName', 'Arial', 'FontSize', 6)
        title(sprintf('frame %d', i))
        pause(0.001)
        %%

%         figure(2); 
%         clf
%         imagesc(currentimg');axis image;colormap('turbo')
%         pause(0.001)
%         xlim([48 90])
%         ylim([68 120])
%         xlabel('X position')
%         ylabel('Y position')
%         xticks(50:10:90)
%         yticks(70:10:120)
%         set(gcf, 'Position', [345,649,144,180])
%         set(gca, 'FontName', 'Arial', 'FontSize', 6)
    end
    J{i} = regions;
    % % Write this frame out to a new video file.
    % writeVideo(outputVideo, getframe(gcf));
end
% % Close video writer
% close(outputVideo);
end

