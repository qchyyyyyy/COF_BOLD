clear;
addpath(genpath('functions'));
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';

[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere(1));
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere(2));

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

%% figure 7c
Condition = 'WM';
PE = 'LR';
load(['data/results/decoding_results/' hemisphere '_decodingAcc_Pattern_individual_' Condition '_' PE '.mat']);
load(['data/results/decoding_results/' hemisphere '_decodingAcc_BOLD_individual_' Condition '_' PE '.mat']);
load(['data/results/decoding_results/' hemisphere '_decodingAcc_dyFC_individual_' Condition '_' PE '.mat']);

Delay = 0:18;
figure;
plot(0:18,mean(decodingAcc_all,2));
p1 = plot(Delay,mean(decodingAcc_all,2),'LineWidth',2);hold on
p2 = plot(Delay,mean(decodingAcc_div,2),'LineWidth',2);hold on
p3 = plot(Delay,mean(decodingAcc_curl,2),'LineWidth',2);hold on
p4 = plot(Delay,mean(decodingAcc_dyFC,2),'LineWidth',2);hold on
p5 = plot(Delay,mean(decodingAcc_BOLD,2),'LineWidth',2);hold on
p6 = plot(Delay,repmat(0.3572,1,19),'--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1.5);hold on
ll = legend([p1,p2,p3,p4,p5,p6],'all patterns','div patterns','curl patterns','dynamic FC','BOLD signals','chance level');
set(ll,'AutoUpdate','off');

pp = patch([Delay flip(Delay)],[mean(decodingAcc_all,2)'+std(decodingAcc_all,[],2)' flip(mean(decodingAcc_all,2)-std(decodingAcc_all,[],2))'],...
    p1.Color,'edgecolor','none');hold on
alpha(pp,0.05);
pp = patch([Delay flip(Delay)],[mean(decodingAcc_div,2)'+std(decodingAcc_div,[],2)' flip(mean(decodingAcc_div,2)-std(decodingAcc_div,[],2))'],...
    p2.Color,'edgecolor','none');hold on
alpha(pp,0.05);
pp = patch([Delay flip(Delay)],[mean(decodingAcc_curl,2)'+std(decodingAcc_curl,[],2)' flip(mean(decodingAcc_curl,2)-std(decodingAcc_curl,[],2))'],...
    p3.Color,'edgecolor','none');hold on
alpha(pp,0.05);
pp = patch([Delay flip(Delay)],[mean(decodingAcc_dyFC,2)'+std(decodingAcc_dyFC,[],2)' flip(mean(decodingAcc_dyFC,2)-std(decodingAcc_dyFC,[],2))'],...
    p4.Color,'edgecolor','none');hold on
alpha(pp,0.05);
pp = patch([Delay flip(Delay)],[mean(decodingAcc_BOLD,2)'+std(decodingAcc_BOLD,[],2)' flip(mean(decodingAcc_BOLD,2)-std(decodingAcc_BOLD,[],2))'],...
    p5.Color,'edgecolor','none');hold on
alpha(pp,0.05);


xlim([0 18]);
ylim([0.2 1]);

xlabel('Delay (TR)','FontSize',14);
yy = ylabel(gca,'Decoding accuracy','FontSize',16);
% set(yy,'position',[-29.304508499630458,28.60352422907491,1.4e-14])
set(gca,'LineWidth',1);

set(gcf,"Position",[341.5,-196,873,304]);