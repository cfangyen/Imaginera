clear;
cd('D:\C elegans Team Dropbox\HONGFEI JI\Dian Hongfei shared\worm analysis\Manuscript\Individual Figures')
load('Sample_HuMSER_longtrace1.mat')
TBends = [0 6 16 23 61 71 82 90 133 158 175 188 212 224 234 249 264 307 318,...
          330 340 354 368 379 389 398 409 419 430 439 450 457 469 487 495,...
          502 510 517 525 533 541 550 558 565 574 582 592 599 607 613 620,...
          628 635 642 648 655 662 669 676 684 691 700 706 715 722 731 739,...
          746 751];
TBends = (TBends-1)/30;
NumT = length(F_seq);
T = (0:(NumT-1))/30;
figure(1); clf;
hold on
for i = 1 : numel(TBends)-1
    currx = [TBends(i) TBends(i+1) TBends(i+1) TBends(i)];
    curry = [-.1 -.1 .1 .1];
    if mod(i, 2) == 1
        patch(currx, curry, 'red', 'EdgeColor', 'none')
    else
        patch(currx, curry, 'blue', 'EdgeColor', 'none')
    end
end
plot(T, F_seq)
hold off
xlim([T(1) T(end)])
ylim([-.1 .1])
yticks([-.1 0 .1])
yticklabels({'-0.1', '0.0', '0.1'})
xticks(15:1:25)
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'})
xlabel('Time (s)')
ylabel(sprintf('Covariance value'))
% yticks([])
xlim([15 25])
ax = gca;
ax.YAxis(1).Exponent = 0;
set(gcf, 'Position', [46,737,380,120])
set(gca, 'FontName', 'Arial', 'FontSize', 6, 'Box', 'off')