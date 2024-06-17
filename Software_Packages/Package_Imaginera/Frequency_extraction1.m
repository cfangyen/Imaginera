function [F0, F0_L, F0_U, F_seq] = Frequency_extraction1(C, fps, delta, movmean_windowratio)
%FREQUENCY_EXTRACTION
% % Compute frequency using peakfinding
N = size(C, 1);
autoC_norm = zeros(1, N);

for i = 1 : N
    autoC_norm(i) = mean(diag(C, i-1), 'omitnan');
end

V = autoC_norm;
szwin = ceil(fps*movmean_windowratio);
Vfil = smoothdata(V, 'movmean', szwin);
V = Vfil;
%%%
if sum(isfinite(V))>=6
    [~, imax0, ~, pmax0] = findpeaks(V);
    [~, imin0, ~, pmin0] = findpeaks(-V);
    imax0   = cat(2, 1, imax0);
    pmax0   = cat(2, 2, pmax0);
    % Identify prominent peaks and troughs
    idx_highpromax = pmax0>=delta;
    idx_highpromin = pmin0>=delta;
    % Identify connected peak/trough chains
    bwconn_imax = bwconncomp(idx_highpromax);
    bwconn_imin = bwconncomp(idx_highpromin);
    imaxidx_list = bwconn_imax.PixelIdxList;
    iminidx_list = bwconn_imin.PixelIdxList;
    imax_list = cell(size(imaxidx_list));
    imin_list = cell(size(iminidx_list));
    for i = 1 : numel(imaxidx_list)
        currimaxidx = imaxidx_list{i};
        imax_list{i} = imax0(currimaxidx);
    end
    for i = 1 : numel(iminidx_list)
        curriminidx = iminidx_list{i};
        imin_list{i} = imin0(curriminidx);
    end
    % Identify peak-trough couples
    ptcouples = {};
    kk = 1;
    for i = 1 : numel(imax_list)
        currimax_chain = imax_list{i};
        currimax_range = currimax_chain(1):currimax_chain(end);
        if ~isempty(imin_list)
            for j = kk : numel(imin_list)
                currimin_chain = imin_list{j};
                currimin_range = currimin_chain(1):currimin_chain(end);
                currintersect = intersect(currimax_range, currimin_range);
                if ~isempty(currintersect)
                    currcouple = {currimax_chain, currimin_chain};
                    kk = kk+1;
                    break;
                else
                    currcouple = [];
                end
            end
        else
            break;
        end
        ptcouples = cat(1, ptcouples, currcouple);
    end
    %%%
    T0_com = [];
    for ii = 1 : size(ptcouples, 1)
        currimax_chain = ptcouples{ii, 1};
        currimin_chain = ptcouples{ii, 2};
        currpt_chain = sort(cat(2, currimax_chain, currimin_chain));
        T0_com = cat(2, T0_com, diff(currpt_chain)*2);
    end

    if numel(T0_com)>=2
        T0 = mean(T0_com)/fps;
        T0_L = prctile(T0_com, 25)/fps;
        T0_U = prctile(T0_com, 75)/fps;
        F0 = 1/T0;
        F0_L = 1/T0_U;
        F0_U = 1/T0_L;
    else
        F0 = NaN;
        F0_L = NaN;
        F0_U = NaN;
    end
else
    F0   = NaN;
    F0_L = NaN;
    F0_U = NaN;
    % ipk = NaN;
end
F_seq = autoC_norm';

end


