% Worm shape analysis. This version is modified for analyzing data obtained
% from worm locomotion on WormTunnel microfluidic chambers.
%
% First, periods during which a worm is undulating without constraint will be
% collected and analyzed to generate a limit cycle for normal undulatory
% dynamics.
%
% Next, periods during which a worm is moving under constraint (usually at
% the middle of body) will be collected and analyzed to generate phase
% dynamics during those periods
%
% Finally, combining the analyzed results from the above two steps, a
% generalized compensatory factor kymogram will be generated in a 2D
% heatmap form, which represents a function                                                                                                                                                                                                         ion of time and body coordinate.
% Also, information of the body portion that is being constrained will be
% reflected on the kymogram.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% START MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
clc; close all; clear

global start_illum end_illum prefix pathname filename
global conc pix_per_mm wormthreshold isie decim filsize spline_p initials


% SScale = [0.03, 0.1, 0.3, 1];
% TScale = [1/15 1/10 1/6 1/3 1/2 1];

resizefactor = 0.4;
tscale = 1;
wormlabel  = 1;
pix_per_mm = 198.83;
prefix;
pathname;

curvlim = 0.2;
domovie = 0;
issavefiles = 0;
iscontinue = 1;

option1 = 'Yes (AVI only)';
option2 = 'Yes (MAT)';

fname  = questdlg('Pre or Post-const undulation?','','Pre-','Post-','All-', 'All-');
button = length(questdlg('Load new data?','',option1,option2,'No', option1));
CURV_all = {};
dCURV_all = {};
LENG_all = {};
AREA_all = {};
pathname_all = {};
if button == length(option2)
    disp('options2');
    [filename,pathname] = uigetfile({'*.mat'});
    matfname  = [fname filename(1:6)];
    load([pathname filename]);
elseif button == length(option1)
    disp('option1');
    do_dialog = 1;
    %% Identify and analyze control periods
    while iscontinue == 1
        
    if do_dialog
        try
            cd(pathname);
        catch
            pathname = pwd;
        end
        [filename,pathname]  = uigetfile('*.avi', 'Select File');
        pathname_all = [pathname_all, pathname];
        matfname  = [fname filename(1:6)];
        vidObj = VideoReader(fullfile(pathname,filename));
        NumFrames = vidObj.NumFrames;
        fps       = vidObj.FrameRate;
% fps       = 23;
        if isempty(conc)
            conc = 0;
        end
        if isempty(wormlabel)
            wormlabel = 1;
        end
        if isempty(pix_per_mm)
            pix_per_mm = 1;
        end
        if isempty(wormthreshold)
            wormthreshold = 0.10;
        end
        if isempty(isie)
            isie = [1, NumFrames];
        end
        if isempty(decim)
            decim = 1;
        end
        if isempty(filsize)
            filsize = 0.2;
        end
        if isempty(start_illum)
            start_illum = 1;
        end
        if isempty(end_illum)
            end_illum = 1;
        end
        if isempty(spline_p)
            spline_p = 0.01;
        end
        if isempty(initials)
            initials = {'JHF'};
        end
        
        fields={'conc','wormlabel', 'fps','pixels per mm','Worm image threshold',...
            'istart/iend (use";"if multiple)', 'Decimation (1=none))',...
            'Filter size / diameter', 'start_illum','end_illum', ...
            'spline fit parameter', 'Make movie?', 'Your initials', 'Save files?'};
        if exist('isie', 'var')
            answer = inputdlg(fields, 'Cancel to clear previous', 1, ...
                {num2str(conc),num2str(wormlabel),num2str(fps),num2str(pix_per_mm),num2str(wormthreshold),...
                mat2str(isie), num2str(decim),...
                num2str(filsize),num2str(start_illum),num2str(end_illum), ...
                num2str(spline_p), num2str(domovie), initials{1}, num2str(issavefiles)});
        else
            answer = inputdlg(fields, '', 1);
        end
        
        if isempty(answer)
            pause;
        end
        
        conc = str2double(answer{1});
        wormlabel = str2double(answer{2});
        fps = str2double(answer{3});
        pix_per_mm = str2double(answer{4});
        wormthreshold = str2double(answer{5});
        isie = Str2Mat(answer{6});
        decim = str2double(answer{7});
        filsize = str2double(answer{8});
        start_illum = str2double(answer{9});
        end_illum = str2double(answer{10});
        spline_p = str2double(answer{11});
        domovie = str2double(answer{12});
        initials = answer(13);
        issavefiles = str2double(answer{14});
    end
    fileID = fopen(fullfile(pathname, 'timestemp for ctrl.txt'),'a+');
    fprintf(fileID, [filename,': ', answer{6},'\n']);
    fclose(fileID);
    nperiods   = size(isie, 1);
    curv_ctrl  = cell(nperiods, 1);
    dcurvdt_ctrl = cell(nperiods, 1);
    angle_ctrl = cell(nperiods, 1);
    len_ctrl   = cell(nperiods, 1);
    area_ctrl  = cell(nperiods, 1);
    %  Compute undulatory variables for control groups
    tic
    for kk = 1 : nperiods
        do_const = 0;
        thisperiod = isie(kk, :);
        options = {conc, wormlabel, fps, pix_per_mm, wormthreshold,...
                   thisperiod, decim, filsize, start_illum, end_illum,...
                   spline_p, domovie, initials, pathname, filename, do_const, issavefiles};
        [curv_ctrl{kk},dcurvdt_ctrl{kk}, angle_ctrl{kk}, len_ctrl{kk}, area_ctrl{kk}]...
            = WORMSHAPE_MAINCALCULATION(vidObj, options, resizefactor, tscale);
    end
    toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(2); clf
%     KK = curv_ctrl{1};
%     NumT = size(KK, 1); T = (0:(NumT-1))/fps;
%     NumS = size(KK, 2); S = linspace(0,1, NumS);
%     imagesc(T,S, KK')
%     colormap(cmap_redblue(0.7))
%     xlabel('Time (s)')
%     ylabel(sprintf('Body coordinate\n(head = 0; tail = 1)'))
% %     clim([-8 8]);
% %     colorbar
%     yticks([0 .5 1])
%     yticklabels({'0.0', '0.5', '1.0'})
%     set(gcf, 'Position', [295,754,300,120])
%     set(gca, 'FontName', 'Arial', 'FontSize', 6)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % flip some periods where head/tail misidentified
%     KK = curv_ctrl{1};
    for i = 1 : nperiods
        K = curv_ctrl{i};
        dKdt = dcurvdt_ctrl{i};
        ang  = angle_ctrl{i};
        % debug plot
        figure(10); clf
        imagesc(K)
        colormap(cmap_redblue(0.7))
        caxis([-25 25])
        colorbar
        hold on
        %%%%%%%%%%%%%%%%%
%         answer = length(questdlg('Need to flip some period?', '','Yes', 'No', 'No'));
        %%%%%%%%%%%%%%%%%
%         if answer == 3
%             title('Indicate the period that need to be flipped')
%             % flip curvature and dKdt
%             [~, flpy1] = ginput(1);
%             flpy1 = max([floor(flpy1) 1]);
%             line([1 100], [flpy1 flpy1],'Color','white','LineStyle','--')
%             [~, flpy2] = ginput(1);
%             flpy2 = min([floor(flpy2) size(K,1)]);
%             line([1 100], [flpy2 flpy2],'Color','white','LineStyle','--')
%             K2flip = K(flpy1 : flpy2, :);
%             dKdt2flip = dKdt(flpy1 : flpy2, :);
%             ang2flip  = ang(flpy1 : flpy2, :);
%             K(flpy1 : flpy2, :) = flip(K2flip,2);
%             dKdt(flpy1 : flpy2, :) = flip(dKdt2flip,2);
%             ang(flpy1 : flpy2, :) = flip(ang2flip,2);
%         end
        %%%%%%%%%%%%%%%%%
        % updata data
        curv_ctrl{i}  = K;
        dcurvdt_ctrl{i} = dKdt;
        angle_ctrl{i}  = ang;
        % show updated plot
        figure(10); clf
        imagesc(K)
        colormap(cmap_redblue(0.7))
        colorbar
        pause(1)
    end
     % end fliping loop
    CURV_all = [CURV_all; curv_ctrl];
    dCURV_all = [dCURV_all; dcurvdt_ctrl];
    LENG_all = [LENG_all; len_ctrl];
    AREA_all = [AREA_all; area_ctrl];
    Cbutton = length(questdlg('Continue?','','Yes','No','Yes'));
    if Cbutton == 3
        iscontinue = 1;
    else
        iscontinue = 0;
    end
    close all
    end
end
%% Calculate the average normal phase plot from control groups
% fps = 30;
fps_current = fps * tscale;
[Kc_all, dKdtc_all, Kp_data, dKdtp_data, T0_avg]...
    = Generalized_cfactor_for_microfluidics(CURV_all, [], fps_current);


numsamplepts = 100;
numcurvpts   = 100;
a  = .15; c = a * T0_avg;
Zc = Kc_all + 1i*c*dKdtc_all;
Pc = unwrap(angle(Zc), [], 2);
[~, Sc] = meshgrid(1:numsamplepts, 1:numcurvpts);
% Generate interpolant (in a bulk manner)
FR = scatteredInterpolant(Pc(:), Sc(:), Zc(:), 'linear', 'nearest');

%%  Plotting data
curvrgn_analysis = 1 : numcurvpts;
% 3-D phase portrait plot of normal undulation
figure(1); clf
TX_mesh = meshgrid(1 : numsamplepts, curvrgn_analysis);
hold on
plot3(Kc_all', dKdtc_all', TX_mesh)  % isosegmental line
plot3(Kc_all,  dKdtc_all,  TX_mesh') % isophasic line
hold off
view([30,30])
xl = xlim;
yl = ylim;
zlim([5 95])
xlabel('K')
ylabel('dKdt')
zlabel('Body coordinate')
set(gca, 'FontSize', 12, 'Position', [0.13,0.11,0.775,0.815])

% GUI with interactive response-plot updates for phase portrait plots
s  = 30;
f  = figure(2); clf
ax = axes('Parent',f,'position',[0.2 0.25  0.65 0.65], 'PlotBoxAspectRatio', [1,0.81,0.75]);
faseplot2  = @(ax, s) phasePlot2(curvrgn_analysis, Kc_all', dKdtc_all', s, ax, xl, yl);
faseplot2(ax, s);
b = uicontrol('Parent',f,'Style','slider','Position',[81,34,419,23],...
              'value',s, 'min',5, 'max',95);
bgcolor = f.Color;
uicontrol('Parent',f,'Style','text','FontSize',12,'Position',[45,37,30,20],...
                'String','Head','BackgroundColor',bgcolor);
uicontrol('Parent',f,'Style','text','FontSize',12,'Position',[505,37,30,20],...
                'String','Tail','BackgroundColor',bgcolor);
uicontrol('Parent',f,'Style','text','FontSize',12,'Position',[240,10,100,23],...
                'String','Body coordinate','BackgroundColor',bgcolor);
            
b.Callback = @(es,ed) faseplot2(ax, es.Value);
%% save the NormalUndulation.mat file to all visited folders

% fname4save = fullfile(pathname_all{1}, [matfname sprintf('imm Normal_S%.2fT%.2f.mat', resizefactor, tscale)]);
% save(fname4save, 'Kc_all', 'dKdtc_all', 'Zc', 'FR', 'CURV_all', 'fps', 'T0_avg', 'len_ctrl')


%%
