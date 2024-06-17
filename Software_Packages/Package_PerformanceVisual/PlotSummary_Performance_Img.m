% Define data
foldername = uigetdir('Select directory with .mat files to process');
dirname = foldername;
load(fullfile(dirname, 'Robustness_test2.mat'), 'R2H', 'R2H2', 'R2', 'Time_cost', 'R3')
CC = linspecer(4);
%%
T_res = [2 3 5 10 15 30]./1.7;
figure(1); clf;
R2H3 = flip(R2H2, 2);
hold on
plot(T_res, R2H3(4,:), 'Color', CC(4,:), 'Marker','diamond')
plot(T_res, R2H3(3,:), 'Color', CC(3,:), 'Marker','square')
plot(T_res, R2H3(2,:), 'Color', CC(2,:), 'Marker','*')
plot(T_res, R2H3(1,:), 'Color', CC(1,:), 'Marker','.')
hold off
xticks([.5 T_res 20])
xticklabels({'', '1.17', '1.76', '2.94', '5.88', '8.82', '17.65', ''})
yticks(0:.2:1)
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xlim([1 20])
% ylim([0 1])
ylabel('R-squared value')
xlabel('Frames per moving cycle')
legend({'150', '20', '6', '4'})
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'XScale', 'log')
set(gcf, 'Position', [295,754,255,150])
%% Accuracy among Trad HuMSER VidCov (varied temporal resolution under 20 px/L spatial resolution)
figure(2); clf;
bar(flip(R2, 2)', 'EdgeColor','none');
yticks([0 .2 .4 .6 .8 1])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ylabel('R-squared value')
xlabel('Frames per moving cycle')
legend({'Imaginera', 'Worm segmentation', 'Video covariance'})
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off')
set(gcf, 'Position', [295,754,200,150])
%% Accuracy among Trad HuMSER VidCov (varied temporal resolution under 6 px/L spatial resolution)
figure(3); clf;
bar(flip(R3, 2)', 'EdgeColor','none');
yticks([0 .2 .4 .6 .8 1])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ylabel('R-squared value')
xlabel('Frames per moving cycle')
legend({'Imaginera', 'Worm segmentation', 'Video covariance'})
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off')
set(gcf, 'Position', [295,754,200,150])
%% Comparing computational complexity
figure(4); clf;
bar(Time_cost, 'EdgeColor','none');
% yticks([0 .2 .4 .6 .8 1])
% yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
% ylim([.01 1])
% yticks(10:100:500)
ylabel('Time to analyze 10k frames (sec)')
xlabel('Frames per moving cycle')
set(gca, 'XGrid', 'off', 'YGrid', 'off', 'FontName', 'Arial','FontSize', 6, 'Box', 'off','YScale','linear')
set(gcf, 'Position', [295,754,150,150])