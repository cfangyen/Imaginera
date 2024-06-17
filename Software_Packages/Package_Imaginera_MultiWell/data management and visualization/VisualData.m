sss = 1;
TT = struct([]);
%%
path_parent = 'C:\Users\fffei\C elegans Team Dropbox\HONGFEI JI\Paper\Egg-laying\Data\Burst MultWell';
path_expt = '2023-09-11 N2-MT14984_5-HT_6cdn';
group1 = 'wt';
group2 = 'tph-1';

cd(fullfile(path_parent, path_expt));
matdir = dir("*_DATA.mat");
load(matdir(1).name)

% Histogram plot
m = 100;
cm_parula=fake_parula(m);
cm_magma=magma(m);
cm_inferno=inferno(m);
cm_plasma=plasma(m);
cm_viridis=viridis(m);
NumT = size(S,1)/8;

S = S(1:NumT, 1);

%% Frequency Activity
FrqAll = cat(3, S.Freq);
PrdAll = cat(3, S.Periodicity);
Frq1_max = FrqAll(1:4,1:4,:);  Frq1_max = reshape(Frq1_max, [16 NumT]);
Frq1_haf = FrqAll(1:4,5:8,:);  Frq1_haf = reshape(Frq1_haf, [16 NumT]);
Frq1_min = FrqAll(1:4,9:12,:); Frq1_min = reshape(Frq1_min, [16 NumT]);
Frq2_max = FrqAll(5:8,1:4,:);  Frq2_max = reshape(Frq2_max, [16 NumT]);
Frq2_haf = FrqAll(5:8,5:8,:);  Frq2_haf = reshape(Frq2_haf, [16 NumT]);
Frq2_min = FrqAll(5:8,9:12,:); Frq2_min = reshape(Frq2_min, [16 NumT]);
Prd1_max = PrdAll(1:4,1:4,:);  Prd1_max = reshape(Prd1_max, [16 NumT]);  Frq1_max = Frq1_max.*(Prd1_max>.3);
Prd1_haf = PrdAll(1:4,5:8,:);  Prd1_haf = reshape(Prd1_haf, [16 NumT]);  Frq1_haf = Frq1_haf.*(Prd1_haf>.3);
Prd1_min = PrdAll(1:4,9:12,:); Prd1_min = reshape(Prd1_min, [16 NumT]);  Frq1_min = Frq1_min.*(Prd1_min>.3);
Prd2_max = PrdAll(5:8,1:4,:);  Prd2_max = reshape(Prd2_max, [16 NumT]);  Frq2_max = Frq2_max.*(Prd2_max>.3);
Prd2_haf = PrdAll(5:8,5:8,:);  Prd2_haf = reshape(Prd2_haf, [16 NumT]);  Frq2_haf = Frq2_haf.*(Prd2_haf>.3);
Prd2_min = PrdAll(5:8,9:12,:); Prd2_min = reshape(Prd2_min, [16 NumT]);  Frq2_min = Frq2_min.*(Prd2_min>.3);

FrqAll = cat(1, Frq1_min, Frq1_haf, Frq1_max,...
                Frq2_min, Frq2_haf, Frq2_max);

Frq1_avg = [mean(Frq1_min,'all','omitnan'), mean(Frq1_haf,'all','omitnan'), mean(Frq1_max,'all','omitnan')];

Frq2_avg = [mean(Frq2_min,'all','omitnan'), mean(Frq2_haf,'all','omitnan'), mean(Frq2_max,'all','omitnan')];

Frq1_sem = [std(Frq1_min,0,'all','omitnan'), std(Frq1_haf,0,'all','omitnan'), std(Frq1_max,0,'all','omitnan')]./sqrt(16*NumT);

Frq2_sem = [std(Frq2_min,0,'all','omitnan'), std(Frq2_haf,0,'all','omitnan'), std(Frq2_max,0,'all','omitnan')]./sqrt(16*NumT);

T = (0:(NumT-1))/30;
% % Representative traces
CC = linspecer(3);
VFrq1_min = Frq1_min(4,:);
VFrq1_haf = Frq1_haf(4,:);
VFrq1_max = Frq1_max(4,:);
figure(7); clf; hold on
plot(T, VFrq1_min, 'Color', CC(1,:), "LineWidth", 1, 'LineStyle', '-')
plot(T, VFrq1_haf, 'Color', CC(2,:), "LineWidth", 1, 'LineStyle', '--')
plot(T, VFrq1_max, 'Color', CC(3,:), "LineWidth", 1, 'LineStyle', '-.')
hold off
xlabel('Time (hr)')
xticks([0 1 2 3 4 5])
ylabel('Frequency (Hz)')
set(gca, 'FontSize', 6)
leg = legend({'0 mM', '6.5 mM', '13 mM'}, 'Location','northeastoutside');
title(leg,'5-HT conc.')
set(gcf, "Position", [418,434,600,150])

figure(8); clf; hold on
imagesc(T, G, FrqAll)
plot([T(1), T(end)], [16.5 16.5], '--c', 'LineWidth', 1);
plot([T(1), T(end)], [32.5 32.5], '--c', 'LineWidth', 1);
plot([T(1), T(end)], [48.5 48.5], '-c', 'LineWidth', 2);
plot([T(1), T(end)], [64.5 64.5], '--c', 'LineWidth', 1);
plot([T(1), T(end)], [80.5 80.5], '--c', 'LineWidth', 1);
colormap(cm_plasma)
title(sprintf('Frequency %s vs %s', group1, group2))
xlabel('Time (hr)')
xticks([0 1 2 3 4 5])
yyaxis left
ylabel('5HT concentration (mM)')
yticks([8 24 40 56 72 88])
xlim([T(1) T(end)])
ylim([1 96])
yticklabels({sprintf('0'), sprintf('6.5'), sprintf('13'),...
             sprintf('0'), sprintf('6.5'), sprintf('13'), '.'})
set(gca, "YDir", "reverse")
yyaxis right
ylabel('Strain')
yticks([24 72])
ylim([1 96])
yticklabels({group1, group2})
set(gca, "YDir", "reverse")
set(gca, 'FontSize', 12, 'Color', 'k')

% % Barplot
figure(9); clf; 
hold on
Y  = [Frq1_avg(1,:); Frq2_avg(1,:)];
YE = [Frq1_sem(1,:); Frq2_sem(1,:)];
ngroups = size(Y, 1);
nbars = size(Y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
bar([Frq1_avg(1,:); Frq2_avg(1,:)])
for i = 1:nbars
    X = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(X, Y(:,i), YE(:,i), '.');
end
legend({'0 mM 5-HT', '6.5 mM', '13 mM'});
xticks([1 2])
xticklabels({group1 group2})
ylabel('Frequency (Hz)')
set(gca, 'FontSize', 12)