function [Kc_all, dKdtc_all, Kp_data, dKdtp_data, T0_avg]...
    = Generalized_cfactor_for_microfluidics(curv_ctrl, curv_const, fps)
% GENERALIZED_CFACTOR defines and calculates a generalized factor to
% quantify the compensation effects of anterior curvature which is induced
% by (mechanically) perturbing the middle curvature during the forward 
% locomotion of a worm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-09-24-20-%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Modified the algorithm for computing generalized compensatory factor of
%   the normal orbits over body coordinates.
%   Specific steps:
%   1. Define the time window(s) of normal undulations.
%   2. Calculate the average normal phase plot where phase angles are
%   defined using peak-finding method
%   3. I divided the body into 5 different segments with equal length and
%   calculate phase plot for each segments individually, since curvature
%   dynamics varies across body coordinates.
%   Specifically, body is divided into segments 10%, 30%, 50%, 70%, and 90%
%   and phase plots elsewhere are computed through interpolation or
%   extrapolation.

Nc = numel(curv_ctrl);
numcurvpts = size(curv_ctrl{1}, 2);
numsamplepts = 100;
curvrgn_analysis = 1 : numcurvpts;
curvrgn_sample = 10:5:90;
curvrgn_rd = 2;
numsegs = numel(curvrgn_sample);
Kc_avg_all = zeros(numsamplepts, numsegs);
dKdtc_avg_all = zeros(numsamplepts, numsegs);
T_all = [];
for ii = 1 : numsegs
    curvrgn = curvrgn_sample(ii);
    fprintf('Computing normal phase plot at region %.0f%%. ', curvrgn)
    Kc_data = {};
    dKdtc_data = {};
    Ntrials_c  = 0;
    curvrgn_win = curvrgn-curvrgn_rd : curvrgn+curvrgn_rd;
    for i = 1:Nc
        curvdatafiltered = curv_ctrl{i};
        v   = mean(curvdatafiltered(:,curvrgn_win),2);
        % Find all peaks
        if anynan(v)
            continue
        end

        [imax, imin] = C2_get_curvature_peaks(v,1);
        % Find the peaks in the period which is not affected by illumination
        [imax, imin] = verify_extrema(v, imax, imin);
        imax = unique(imax);
        imin = unique(imin);
        imax(v(imax)<=0) = [];
        imin(v(imin)>=0) = [];
        T0 = mean([diff(imax); diff(imin)])/fps;
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        zx = zci(v);
        num_max = numel(imax);
        for j = 1 : num_max - 1
            imaxs = imax(j);
            imaxe = imax(j+1);
            current_zx = zx(zx>imaxs & zx<imaxe);
            if numel(current_zx)~=2
                continue;
            end
            K = v(imaxs: imaxe);
            dKdt = gradient(K)*fps;
            dKdt(1) = 0;
            dKdt(end) = 0;
            Kc_data = [Kc_data; K];
            dKdtc_data = [dKdtc_data; dKdt];
            Ntrials_c = Ntrials_c + 1;
        end
        % record the period into matrix T
        current_seg = curvrgn_sample(ii);
        if current_seg>=20 && current_seg<=40
            T_all   = [T_all, T0];
        end
    end
    numhfsamplepts = numsamplepts/2;
    N_cycs = numel(Kc_data);
    Kc_resc_data = zeros(N_cycs, numsamplepts);
    dKdtc_resc_data = zeros(N_cycs, numsamplepts);
    for i = 1 : N_cycs
        v    = Kc_data{i};
        dvdt = dKdtc_data{i};
        [~, imin] = min(v);
        if length(imin)>1
            tmpdist2midpt = abs(imin - length(v)/2);
            [~,tmpI] = min(tmpdist2midpt);
            imin = imin(tmpI);
        end
        % rescale Ks and dKdts half cycle by half cycle
        vhf1 = v(1:imin); vhf2 = v(imin:end);
        qx1  = numel(vhf1) / (numel(vhf1)-1);
        qx2  = numel(vhf2) / (numel(vhf2)-1);
        vhf1_resc = interp1((0:numel(vhf1)-1), vhf1, (numel(vhf1)-1)*(0:numhfsamplepts-1)/(numhfsamplepts-1),'linear');
        vhf2_resc = interp1((0:numel(vhf2)-1), vhf2, (numel(vhf2)-1)*(0:numhfsamplepts-1)/(numhfsamplepts-1),'linear');
        
        dvdthf1 = dvdt(1:imin); dvdthf2 = dvdt(imin:end);
        qy1  = numel(dvdthf1) / (numel(dvdthf1)-1);
        qy2  = numel(dvdthf2) / (numel(dvdthf2)-1);
        dvdthf1_resc = interp1((0:numel(dvdthf1)-1), dvdthf1, (numel(dvdthf1)-1)*(0:numhfsamplepts-1)/(numhfsamplepts-1),'linear');
        dvdthf2_resc = interp1((0:numel(dvdthf2)-1), dvdthf2, (numel(dvdthf2)-1)*(0:numhfsamplepts-1)/(numhfsamplepts-1),'linear');
        
        % combine halfs to get full cycles
        v_rescaled = [vhf1_resc vhf2_resc];
        dvdt_rescaled = [dvdthf1_resc dvdthf2_resc];
        Kc_resc_data(i,:) = v_rescaled;
        dKdtc_resc_data(i,:) = dvdt_rescaled;
    end
    % Calculating averages for current segment
    Kc_resc_avg    = mean(Kc_resc_data,1);
    dKdtc_resc_avg = mean(dKdtc_resc_data,1);
    Kc_avg_all(:, ii) = Kc_resc_avg;
    dKdtc_avg_all(:, ii) = dKdtc_resc_avg;
    fprintf('\n'); % To go to a new line after reaching 100% progress
end
T0_avg = mean(T_all, 'omitnan');
% Adjust the averaged curvature cycle by fixing the minimum point and 
% scaling the near-end points so that it will equal to the negative 
% amplitude
for ii = 1 : numsegs
    v  = Kc_avg_all(:,ii);
    [~,pf] =  min(v); pf = pf(1);
    vf = abs(v(pf));
    v(1 : pf-1) = v(1 : pf-1) + (pf - (1 : pf-1)')/(pf - 1) * (vf - v(1));
    v(pf : end) = v(pf : end) + ((pf : numel(v))' - pf)/(numel(v) - pf) * (vf - v(end));
    Kc_avg_all(:,ii) = v;
end
% Using interpolation and extrapolation to predict the phase plots on other
% segments of a worm
Kc_all    = interp1(curvrgn_sample, Kc_avg_all',curvrgn_analysis, 'makima')';
dKdtc_all  = interp1(curvrgn_sample, dKdtc_avg_all',curvrgn_analysis, 'makima')';
Kc_all    = Kc_all';
dKdtc_all = dKdtc_all';
%%
%%%%%%%%%%%%%%%%%%%%
if isempty(curv_const)
    Kp_data = [];
    dKdtp_data = [];
else
    % analysis of the contrained group
    Np = numel(curv_const);
    Ntrials_p = 0;
    Kp_data    = {};
    dKdtp_data = {};
    for i = 1:Np
        fprintf('Analyzing trial %d',i)
        curvdatafiltered = curv_const{i};
        
        Ntrials_p = Ntrials_p + 1;
        
        % Analyzing bulk curvature
        Kb    = curvdatafiltered(:,curvrgn_analysis);
        dKdtb = gradient(Kb')' * fps;
        
        Kp_data    = [Kp_data; Kb'];
        dKdtp_data = [dKdtp_data; dKdtb'];
        
        fprintf('\n')
    end
end
end

    function [imax,imin] = verify_extrema(v,imax,imin)
        
        % get the mean amplitudes
        vmax = mean(v(imax));
        vmin = mean(v(imin));
        
        if vmin > vmax
            itemp = imax;
            imax  = imin;
            imin  = itemp;
        end
        
        
    end