function [J, J_noiseremoved] = MSER_feature_extraction_mw(J_feature, fps, nbins, movmean_windowratio)
%MSER_FEATURE_EXTRACTION: Extracting features from MSER objects
numFrames = numel(J_feature);
numFeatures = nbins*2;
J = nan(numFrames, numFeatures);
J_noiseremoved = nan(size(J));
% Tm = nan(numFrames, 1);
RangeOfArea = nan(numFrames, 2); % [min, max]
for i = 1 : numFrames
    currJ = J_feature{i};
    if currJ.Count>=2
        currPixelList = currJ.PixelList;
        minArea = size(currPixelList{1}, 1);
        maxArea = size(currPixelList{end}, 1);
        RangeOfArea(i,:) = [minArea, maxArea];
    end
end
minAreas = RangeOfArea(:, 1);
maxAreas = RangeOfArea(:, 2);
medMinArea = median(minAreas, 'omitnan');
medMaxArea = median(maxAreas, 'omitnan');
X1 = linspace(medMinArea, medMaxArea, nbins);
X1 = unique(X1, "stable");
if numel(X1) == nbins
    for i = 1 : numFrames
        currJ = J_feature{i};
        if currJ.Count>=2
            currNumObj = currJ.length;
            X0 = nan(1, currNumObj);
            E0 = nan(1, currNumObj);
            Or0 = nan(1, currNumObj);
            currOrientation = currJ.Orientation;
            for j = 1 : currNumObj
                % Exclude objects that are too round or too eccentric
                currAxes = currJ(j).Axes;
                a = currAxes(1);
                b = currAxes(2);
                Eccentricity = sqrt(1 - b^2/a^2);
                if Eccentricity<0.05
                    E0(j) = nan;
                    X0(j) = nan;
                    Or0(j) = nan;
                else
                    E0(j) = Eccentricity;
                    X0(j) = size(currJ(j).PixelList,1);
                    Or0(j) = currOrientation(j);
                end
            end
            TanOr0 = tan(Or0);
            T0_clip = atan(TanOr0);
            T0_mean = angle(sum(exp(1i*T0_clip), 'omitnan'));
            T0 = T0_clip - T0_mean;
            [X0C, ia] = unique(X0, 'stable');
            E0C = E0(ia);
            T0C = T0(ia);
            idx2rm = isnan(X0C);
            X0C = X0C(~idx2rm);
            E0C = E0C(~idx2rm);
            T0C = T0C(~idx2rm);
            if numel(E0C)>=2
                E1 = interp1(X0C, E0C, X1, 'linear','extrap');
                T1 = interp1(X0C, T0C, X1, 'linear','extrap');
                J(i, 1:nbins) = E1.*T0_mean;
                J(i, (1:nbins)+nbins) = T1;
%                 Tm(i) = T1(end) - T1(1);
            end
        end
    end
    
    for k = 1 : numFeatures
        V = J(:, k);
        V_filtered = ReduceSignalNoise_mw(V', fps, movmean_windowratio);
        J_noiseremoved(:, k) = V_filtered;
    end
end

end

