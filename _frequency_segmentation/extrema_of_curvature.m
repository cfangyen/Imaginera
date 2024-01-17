% extrema_of_curvature.m
% Anthony Fouad
% Fang-Yen Lab
% 2/26/2015
%
% Calculate the locations where the curvature over time vector v reversed
% direction. Better alternative to extrema.m which is very noisy. 
%
% INPUTS:
%       v           - the curvature over time vector
%
%       thresh      - the curvature value below which data is likely noise. 
%                     0.01 by default. 
%
%       nondupflag  - if 1


function [vmax,imax,vmin,imin] = extrema_of_curvature(v,thresh)

%error('DEPRECATED: Better results obtained using extrema_robust_v2017a.m, which makes use of matlabs built-in findpeaks function');

% Handle inputs
    if nargin < 2; thresh = 0.01; end
    if nargin < 3; nondupflag=0;  end

% Filter out noise
    v(abs(v)<thresh) = 0;

% Break up the signal into positive and negative regions;
    ipos = v>0;
    ineg = v<0;
    
% Number each region
    lpos = bwlabel(ipos);
    lneg = bwlabel(ineg);
    
% Find MAXIMA: the highest value in each positive region
    imax = zeros([max(lpos),1]);
    vmax = zeros([max(lpos),1]);
    
    for i = 1:max(lpos)
        
        % Exactract only this positive region
        vtemp       = v .* (lpos==i);
        
        % Find its maxima
        vmax(i)     = max(vtemp);
        imax(i)     = find(vtemp == max(vtemp),1);
        
    end
    
    
% Find MINIMA: the lowest value in each negative region  
    imin = zeros([max(lneg),1]);
    vmin = zeros([max(lneg),1]);    
  
    for i = 1:max(lneg)
        
        % Exactract only this positive region
        vtemp       = v .* (lneg==i);
        
        % Find its maxima
        vmin(i)     = min(vtemp);
        imin(i)     = find(vtemp == min(vtemp),1);
        
    end    
    
% Clean up results: make sure that there are not two maxs or mins in a row
    
    if nondupflag==1
        % mush all of the i values together
        iall = sort([imin(:) ; imax(:)]);

        % Remember which ones are maximas (+1) or minimas (-1)
        for i = 1:length(iall)
            if v(iall(i))>0; iall(i,2) = 1;
            else  iall(i,2) = -1;
            end
        end

        % Delete duplicate peaks. For each object, take their average location
        d = [diff(iall(:,2));2];
        to_delete = d==0;
        iall(to_delete,:) = [];

        % Now, reassemble the i values into imin and imax
        imin = iall(iall(:,2)==-1,1);
        imax = iall(iall(:,2)==1,1);

        % Finally, recover the actual curvature values at those locations
        vmin = v(imin);
        vmax = v(imax);
    end
% All done!    
end