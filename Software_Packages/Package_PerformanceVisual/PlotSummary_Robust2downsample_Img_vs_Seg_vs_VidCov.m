clc; clear; close all
fs2fo = [2.9, 5.9, 8.8, 17.7];
ParentDir = 'C:\Users\fffei\Dropbox\dropbox2016-2024\Paper\Automatic frequency measurement\worm analysis\Data\Videos_consol';
cd(ParentDir)
SheetNames = 'Robustness to noisy videos.xlsx';
NumCdn = numel(fs2fo);
value = zeros(1, NumCdn);

tempRes = 'tempRes';
Img = 'Img';
Seg = 'Seg';
Cov = 'Cov';

R2_20px = struct(tempRes,fs2fo,Img,value,Seg,value,Cov,value);
R2_06px = struct(tempRes,fs2fo,Img,value,Seg,value,Cov,value);

R2_20px.Img = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'B3:B6');
R2_20px.Seg = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'C3:C6');
R2_20px.Cov = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'D3:D6');

R2_06px.Img = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'E3:E6');
R2_06px.Seg = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'F3:F6');
R2_06px.Cov = readmatrix(SheetNames, 'FileType', 'spreadsheet', 'Sheet', 2,  'Range', 'G3:G6');

R2_Img_20px = R2_20px.Img; R2_Img_20px(R2_Img_20px<0) = 0;
R2_Seg_20px = R2_20px.Seg; R2_Seg_20px(R2_Seg_20px<0) = 0;
R2_Cov_20px = R2_20px.Cov; R2_Cov_20px(R2_Cov_20px<0) = 0;

R2_Img_06px = R2_06px.Img; R2_Img_06px(R2_Img_06px<0) = 0;
R2_Seg_06px = R2_06px.Seg; R2_Seg_06px(R2_Seg_06px<0) = 0;
R2_Cov_06px = R2_06px.Cov; R2_Cov_06px(R2_Cov_06px<0) = 0;

Data_R2_20px = cat(1, R2_Img_20px', R2_Seg_20px', R2_Cov_20px');
Data_R2_06px = cat(1, R2_Img_06px', R2_Seg_06px', R2_Cov_06px');

figure(1); clf;
bar(1:4, Data_R2_20px', 'EdgeColor','none');
yticks([0 .2 .4 .6 .8 1])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xticklabels({'2.9', '5.9', '8.8', '17.7'})
ylabel('R-squared value')
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off')
set(gcf, 'Position', [295,354,580,360])
set(gca, 'FontSize', 15)

figure(2); clf;
bar(1:4, Data_R2_06px', 'EdgeColor','none');
yticks([0 .2 .4 .6 .8 1])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xticklabels({'2.9', '5.9', '8.8', '17.7'})
ylabel('R-squared value')
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off')
set(gcf, 'Position', [295,354,580,360])
set(gca, 'FontSize', 15)
