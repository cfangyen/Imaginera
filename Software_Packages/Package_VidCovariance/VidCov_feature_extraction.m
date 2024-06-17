function [F0, F0_L, F0_U, F_seq] = VidCov_feature_extraction(C, fps, delta)
%VIDCOV_FEATURE_EXTRACTION
N = size(C, 1);
avgC = mean(C, 1, 'omitnan');
stdC = std(C, [], 1, 'omitmissing');
autoC_norm = zeros(1, N);
autoC_orig = zeros(1, N);

for i = 1 : N
    autoC_norm(i) = mean(diag(C, i-1), 'omitnan');
    autoC_orig(i) = mean(diag(C, i-1), 'omitnan')/N*(N-i+1);
end
V = autoC_norm;
szwin = ceil(fps*1/6);
Vfil = smoothdata(V, 'movmean', szwin);
V = Vfil;
NT = numel(V);
%%%
% figure(2); clf
% hold on
% plot(1:NT, V);
% plot(1:NT, Vfil);
% xlabel('Time (s)')
% ylabel('Auto-correlation')
% ylim([-1 1])
%%%
if sum(isfinite(V))>=6
    [pksmax0, imax0, ~, pmax0] = findpeaks(V);
    [pksmin0, imin0, ~, pmin0] = findpeaks(-V);
    pksmax0 = cat(2, 1, pksmax0);
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
    % figure(2); hold on
    % % Indicate all peaks
    % for ii = 1 : numel(pksmax0)
    %     if idx_highpromax(ii) == 1
    %         plot(imax0(ii), pksmax0(ii), 'or')
    %     elseif idx_highpromax(ii) == 0
    %         plot(imax0(ii), pksmax0(ii), 'ok')
    %     end
    % end
    % for ii = 1 : numel(pksmin0)
    %     if idx_highpromin(ii) == 1
    %         plot(imin0(ii), -pksmin0(ii), '<r')
    %     elseif idx_highpromin(ii) == 0
    %         plot(imin0(ii), -pksmin0(ii), '<k')
    %     end
    % end
    % %%% Indicate couples
    % fnh = @sRGB_to_OKLab;
    % rgb = maxdistcolor(size(ptcouples, 1),fnh);
    % for ii = 1 : size(ptcouples, 1)
    %     currimax_chain = ptcouples{ii, 1};
    %     currimin_chain = ptcouples{ii, 2};
    %     plot(currimax_chain, V(currimax_chain), 'Marker','o', 'MarkerFaceColor', rgb(ii,:), 'LineStyle','none', 'MarkerEdgeColor','none')
    %     plot(currimin_chain, V(currimin_chain),'Marker','<',  'MarkerFaceColor', rgb(ii,:), 'LineStyle','none', 'MarkerEdgeColor','none')
    % end
    % set(gcf, 'Position', [235,467,560,420])
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
F_seq{1} = avgC;
F_seq{2} = stdC;
F_seq{3} = autoC_norm'; % This is being used for frequency extraction
F_seq{4} = autoC_orig';
end

