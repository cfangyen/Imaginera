clc; clear; close all
noiselvl  = [0 0.0001, 0.0003, 0.0009, 0.0027, 0.0081, 0.0243 0.0729];
snrValues =   [25.92,  21.62,  17.31,  13.01,  8.71,   4.40,  0.10];
ParentDir = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol';
cd(ParentDir)
SheetNames = 'Robustness to noisy videos.xlsx';
NumCdn = numel(noiselvl);
value = zeros(1, NumCdn);
noiselevel = 'noiselvl';
RMSE  = 'RMSE';
RegR2 = 'RegR2';
RobR2 = 'RobR2';
NumEx = 'NumEx';
Perfm_Img = struct(noiselevel,value,RMSE,value,RegR2,value,RobR2,value,NumEx,value);
Perfm_Seg = struct(noiselevel,value,RMSE,value,RegR2,value,RobR2,value,NumEx,value);
Perfm_Cov = struct(noiselevel,value,RMSE,value,RegR2,value,RobR2,value,NumEx,value);

Perfm_Img.noiselvl = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B2:I2');
Perfm_Img.RMSE     = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B3:I3');
Perfm_Img.RegR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B4:I4');
Perfm_Img.RobR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B5:I5');

Perfm_Seg.noiselvl = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B9:I9');
Perfm_Seg.RMSE     = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B10:I10');
Perfm_Seg.RegR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B11:I11');
Perfm_Seg.RobR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B12:I12');

Perfm_Cov.noiselvl = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B16:I16');
Perfm_Cov.RMSE     = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B17:I17');
Perfm_Cov.RegR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B18:I18');
Perfm_Cov.RobR2    = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 1,  'Range', 'B19:I19');

Perfm_Img_noise = Perfm_Img.RegR2; Perfm_Img_noise(Perfm_Img_noise<0) = 0;
Perfm_Seg_noise = Perfm_Seg.RegR2; Perfm_Seg_noise(Perfm_Seg_noise<0) = 0;
Perfm_Cov_noise = Perfm_Cov.RegR2; Perfm_Cov_noise(Perfm_Cov_noise<0) = 0;

Data_RobR2 = cat(1, Perfm_Img_noise(2:end), Perfm_Seg_noise(2:end), Perfm_Cov_noise(2:end));
figure(1); clf;
bar(snrValues, Data_RobR2', 'EdgeColor','none');
yticks([0 .2 .4 .6 .8 1])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xticks(flip(snrValues))
% xticklabels({'3.3', '7.4', '12.1', '16.8', '21.4', '25.9'})
ylabel('R-squared value')
xlabel('Signal-to-Noise Ratio (dB)')
% set(gca, 'XScale', 'log')
% legend({'Imaginera', 'Segmentation'})
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off')
set(gcf, 'Position', [295,354,800,360])
set(gca, 'FontSize', 15)