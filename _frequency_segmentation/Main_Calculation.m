%% Perform background subtraction for individual videos, save
clear; clc; close all
p = {
    'D:\Dropbox\Paper\Egg-laying\Data\2023-3-18 2hr 2min 15sec 20fps 6 well dop3';...
    };

parentfolder = p{1};
Dir_orgvid = dir(fullfile(parentfolder, 'original', '*.avi'));
wellimgdir = fullfile(parentfolder, 'wellimg');
wellviddir = fullfile(parentfolder, 'wellvid');
mkdir(wellimgdir)
mkdir(wellviddir)
prompt = {'How many wells'};
dlgtitle = '';
dims = [1 20];
definput = {''};
WellNum = inputdlg(prompt,dlgtitle,dims,definput);
WellNum = WellNum{1};
numwells = str2double(WellNum);
numvids  = numel(Dir_orgvid);
separation = 1;
% parpool(12)
for ii = 1 : numvids
    fname = Dir_orgvid(ii).name;
    orgvidfolder = Dir_orgvid(ii).folder;
    tic
    Vid_BkgSubtraction_Output(fname, orgvidfolder, wellimgdir, wellviddir, 'd', numwells, separation);
    toc
end
% delete(gcp('nocreate'))
%% Count eggs from individual bkg images (by calling Python scripts)
% Set the Python executable to the one in the virtual environment
pyenv("Version", "C:\Users\fffei\PycharmProjects\pythonProject\venv\Scripts\python.exe");
% Check the Python version and executable used by MATLAB
current_py_env = pyenv;
disp("MATLAB is using Python version: " + string(current_py_env.Version));
disp("Python executable: " + string(current_py_env.Executable));
% Add the path to the Python module
module_path = fullfile('C:\Users\fffei\Dropbox\Paper\Egg-laying\Code\RCNN\Source_code\mask_rcnn');
if count(py.sys.path, module_path) == 0
    insert(py.sys.path, int32(0), module_path);
end
% Import the Python module
Mask_R_CNN = py.importlib.import_module('Tutorial_Mask_R_CNN_PyTorchOfficial');
Mask_R_CNN = py.importlib.reload(Mask_R_CNN);
% Call the Python module
mfolder_path = fullfile(p{1}, 'wellimg');
Mask_R_CNN.run_main(pyargs('mfolder_path', mfolder_path));
fprintf('Finished')
%% Analyze locomotion activity from pixel difference in burst videos
Dir_mdfvid = dir(fullfile(wellviddir, 'Well*'));
LocoAVG_data = cell(numwells, 1);
LocoSEM_data = cell(numwells, 1);
for i = 1 : numwells
    currentvidpath = fullfile(Dir_mdfvid(i).folder, Dir_mdfvid(i).name);
    matfiledir = dir(fullfile(currentvidpath, '*.mat'));
    numvids  = numel(matfiledir);
    currentwellLocoAVG = zeros(numvids, numel(separation));
    currentwellLocoSEM = zeros(numvids, numel(separation));
    for j = 1 : numvids
        load(fullfile(currentvidpath, matfiledir(j).name), 'LocoAVG', 'LocoSEM')
        currentwellLocoAVG(j,:) = LocoAVG;
        currentwellLocoSEM(j,:) = LocoSEM;
    end
    LocoAVG_data{i} = currentwellLocoAVG;
    LocoSEM_data{i} = currentwellLocoSEM;
end
% % Save data for locomotion dynamics from pixel difference
save(fullfile(parentfolder, 'LocoEgg_data.mat'),...
    'LocoAVG_data', 'LocoSEM_data');
%% Analyze egg-laying activity from Faster R-CNN detection
Dir_wellimg = dir(fullfile(wellimgdir, 'Well*'));
EggAVG_data = cell(numwells, 1);
EggSEM_data = cell(numwells, 1);
for i = 1 : numwells
    fprintf('Analyzing Well %d out of %d\n', i, numwells)
    currentimgpath = fullfile(Dir_wellimg(i).folder, Dir_wellimg(i).name);
    load(fullfile(currentimgpath, 'TOTAL_MASKNUMS.mat'), 'TOTAL_MASKNUMS');
    numimgs = numel(TOTAL_MASKNUMS);
    currentwellEggAVG = zeros(numvids, 1);
    currentwellEggSEM = zeros(numvids, 1);
    numimages_done = 0;
    for j = 1 : numvids
        currentVidImgsdir = dir(fullfile(currentimgpath, sprintf('*fps%04d*', j)));
        numimages_inqueue = numel(currentVidImgsdir);
        currentwellEggAVG(j) = mean(double(TOTAL_MASKNUMS((numimages_done+1) : (numimages_done+numimages_inqueue))));
        currentwellEggSEM(j) = std(double(TOTAL_MASKNUMS((numimages_done+1) : (numimages_done+numimages_inqueue))))/sqrt(numimages_inqueue);
        numimages_done = numimages_done + numimages_inqueue;
    end
    EggAVG_data{i} = currentwellEggAVG;
    EggSEM_data{i} = currentwellEggSEM;
end
% % Save data for locomotion dynamics from pixel difference
save(fullfile(parentfolder, 'LocoEgg_data.mat'),...
    'EggAVG_data', 'EggSEM_data', '-append');
%% Plot egg-laying and locomotion dynamics
load(fullfile(parentfolder, 'LocoEgg_data.mat'))
Genotype = {'dop-3'};
Interval = 2; 
numGroups = numel(Genotype);
% % Plot egg-laying activity
figure(1); clf; hold on
for i = 1 : numwells
    currentWellEggAVG = EggAVG_data{i};
    currentWellEggSEM = EggSEM_data{i};
    numvids = numel(currentWellEggAVG);
    T = ((1 : numvids)-1)*Interval;
    errorbar(T, currentWellEggAVG, currentWellEggSEM, 'LineWidth',2)
end
hold off
xlabel('Time (min)')
ylabel('# of eggs')
set(gca, 'FontSize', 18)
set(gcf, 'Position', [600,500,781,347])
% % Plot locomotion activity
figure(2); clf; hold on
separa2use = 1;
for i = 1 : numwells
    currentwellLocoAVG = LocoAVG_data{i}(:,separa2use);
    currentwellLocoSEM = LocoSEM_data{i}(:,separa2use);
    numvids = numel(currentwellLocoAVG);
    T = ((1 : numvids)-1)*Interval;
    errorbar(T, currentwellLocoAVG, currentwellLocoSEM, 'LineWidth',2)
end
hold off
xlabel('Time (min)')
ylabel('Loco Activity')
set(gca, 'FontSize', 18)
set(gcf, 'Position', [600,500,781,347])
% % Plot time separation vs corresponding locomotion activity
FPS = 20;
SEP = separation.*1/FPS;
LocoAVG_all = cat(3, LocoAVG_data{:});
LocoAVG_norm = LocoAVG_all;
for i = 1 : size(LocoAVG_all,1)
    for j = 1 : size(LocoAVG_all,3)
        LocoAVG_norm(i,:,j) = LocoAVG_norm(i,:,j)./LocoAVG_norm(i,end,j);
    end
end
LocoAVG_wellmerge = mean(LocoAVG_norm, 3);
LocoAVG_merge = mean(LocoAVG_wellmerge, 1);
LocoSEM_merge = std(LocoAVG_wellmerge,1, 1)/sqrt(numvids);
figure(3); clf; hold on
errorbar(SEP, LocoAVG_merge, LocoSEM_merge, 'LineWidth',2)
xlabel('Time separation (sec)')
ylabel('Loco Activity')
set(gca, 'FontSize', 18)
set(gcf, 'Position', [600,500,781,347])
