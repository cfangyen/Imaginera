% C2_get_curvature_peaks.m
% Anthony Fouad
% Fang-Yen Lab
% Summer 2015
%
% Wrapper for extrema_of_curvature.m. Extracts the peaks of a worm's
% curvature trace in input argument v. Performs filtering first to remove
% low frequency fluctuations from consideration, and thresholds peaks below
% 0.05 amplitude (out of +/- 0.5)
%
% INPUTS:
% 
%       v       -   vector of curvature values
%
%       precise -   Whether to return the true peaks (1) or the approximate
%                   peaks from the smoothened curve used for suppressing
%                   noise (0). 
%
%       thresh  -   maxima/minima must be above/below this/-this threshold
%                   to be considered. Note that the sinusoid is rescaled to
%                   lie between [-0.5,0.5]




function [imax,imin]=C2_get_curvature_peaks(v,precise,thresh)

% Handle inputs
    if nargin < 2; precise = 1; end
    if nargin < 3; thresh = 0.1;end

% Get low frequency components of the curve using a cubic fit
    x       = 0:length(v)-1; x = x(:); v = v(:);
    fit     = createSplineFit(x,v,3.353875366298346E-5);
    vlf     = fit(x);
    vhf     = v(:)-vlf(:);
    
% Exclude outlier regions before normalizing
    loval = prctile(vhf,5);
    hival = prctile(vhf,95);
    vhf(vhf>hival) = hival;
    vhf(vhf<loval) = loval;
    
% Normalize the high frequency components
    vhf = normalize21(vhf)-0.5;

% Get the maxima and minima
    [~,imax,~,imin] = extrema_of_curvature(vhf,thresh);

% Refind all maxima if precise is requested
    if precise
        
        % Find extrema in the raw, unsmoothend vector
        [~,imaxraw,~,iminraw] = extrema(v);
        
        % For each extrema calculated robustly, find the nearest raw
        % extrema
        k       = dsearchn(imaxraw,imax);
        imaxnew = imaxraw(k);
        
        k       = dsearchn(iminraw,imin);
        iminnew = iminraw(k);
        
        % Return the precise values
        imax    = imaxnew;
        imin    = iminnew;
    end

end