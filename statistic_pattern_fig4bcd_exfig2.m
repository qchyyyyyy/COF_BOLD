clear;clc
addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks'
[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere{1});
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere{2});

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
%% figure 4b (distribution of dwelling time)
patternNames = {'diverging','converging','CC-spiral','C-spiral'};

load('data/results/dwelltime_results/lh_PatternDwellTime_group_k=2.mat');
original.Dwell_time_Source = Dwell_time_Source;
original.Dwell_time_Sink = Dwell_time_Sink;
original.Dwell_time_CCVortex = Dwell_time_CCVortex;
original.Dwell_time_CVortex = Dwell_time_CVortex;

load('data/results/dwelltime_results/lh_PatternDwellTime_group_k=2_original.mat');
surrogate.Dwell_time_Source = Dwell_time_Source;
surrogate.Dwell_time_Sink = Dwell_time_Sink;
surrogate.Dwell_time_CCVortex = Dwell_time_CCVortex;
surrogate.Dwell_time_CVortex = Dwell_time_CVortex;

%
fontSize = 15;
alphaValue=0.3;

figure;
subplot(2,2,1);
hh = histogram(original.Dwell_time_Source(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.Dwell_time_Source(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title(patternNames{1},'FontSize',fontSize);
xlim([0 30]);ylim([0 0.3]);
ylabel('Fraction');
legend('COF','optical flow');
set(gca,'LineWidth',1);
subplot(2,2,2);
hh = histogram(original.Dwell_time_Sink(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.Dwell_time_Sink(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title(patternNames{2},'FontSize',fontSize);
xlim([0 30]);ylim([0 0.3]);
legend('COF','optical flow');
set(gca,'LineWidth',1);
subplot(2,2,3);
hh = histogram(original.Dwell_time_CCVortex(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.Dwell_time_CCVortex(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title(patternNames{3},'FontSize',fontSize);
xlim([0 30]);ylim([0 0.3]);
ylabel('Fraction');xlabel('Dwelling time of patterns (TR)');
legend('COF','optical flow');
set(gca,'LineWidth',1);

subplot(2,2,4);
hh = histogram(original.Dwell_time_CVortex(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.Dwell_time_CVortex(:),linspace(0,30,31)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title(patternNames{4},'FontSize',fontSize);
xlim([0 30]);ylim([0 0.3]);
ylabel('Fraction');xlabel('Dwelling time of patterns (TR)');
legend('COF','optical flow');
set(gca,'LineWidth',1);

set(gcf,'Position',[458,172,639.5,526.5])
%% distribution of the number of patterns
Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_REST1_LR.mat']);
Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_REST1_LR.mat']);
NumSource = Num_L.NumSource + Num_R.NumSource;
NumSink = Num_L.NumSink + Num_R.NumSink;
NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;

NumSource = NumSource(:,101:1100);
NumSink = NumSink(:,101:1100);
NumCCVortex = NumCCVortex(:,101:1100);
NumCVortex = NumCVortex(:,101:1100);

original.NumSource = NumSource;
original.NumSink = NumSink;
original.NumCCVortex = NumCCVortex;
original.NumCVortex = NumCVortex;

Num_L = load(['data/results/numPattern_results/global/surrogate/' hemisphere{1} '_numPattern_individual_REST1_LR.mat']);
Num_R = load(['data/results/numPattern_results/global/surrogate/' hemisphere{2} '_numPattern_individual_REST1_LR.mat']);
NumSource = Num_L.NumSource + Num_R.NumSource;
NumSink = Num_L.NumSink + Num_R.NumSink;
NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;

NumSource = NumSource(:,101:1100);
NumSink = NumSink(:,101:1100);
NumCCVortex = NumCCVortex(:,101:1100);
NumCVortex = NumCVortex(:,101:1100);

surrogate.NumSource = NumSource;
surrogate.NumSink = NumSink;
surrogate.NumCCVortex = NumCCVortex;
surrogate.NumCVortex = NumCVortex;
%
fontSize = 15;
alphaValue=0.3;
Interval = 100;

figure;
clf
subplot(2,2,1);
hh = histogram(original.NumSource(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.NumSource(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title('Source','FontSize',fontSize);
xlim([0 Interval+0.5]);ylim([0 0.08]);
ylabel('Fraction');
legend('original','surrogate');
set(gca,'LineWidth',1);

subplot(2,2,2);
hh = histogram(original.NumCCVortex(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.NumCCVortex(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title('CC-spiral','FontSize',fontSize);
xlim([0 Interval+0.5]);ylim([0 0.08]);
legend('original','surrogate');
set(gca,'LineWidth',1);

subplot(2,2,3);
hh = histogram(original.NumSink(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.NumSink(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title('Sink','FontSize',fontSize);
xlim([0 Interval+0.5]);ylim([0 0.08]);
ylabel('Fraction');xlabel('The number of patterns per frame');
legend('original','surrogate');
set(gca,'LineWidth',1);

subplot(2,2,4);
hh = histogram(original.NumCVortex(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
hold on
hh = histogram(surrogate.NumCVortex(:),linspace(0,Interval,Interval+1)+0.5,'EdgeAlpha',0,'Normalization','probability');
alpha(hh,alphaValue);
title('C-spiral','FontSize',fontSize);
xlim([0 Interval+0.5]);ylim([0 0.08]);
xlabel('The number of patterns per frame');
legend('original','surrogate');
set(gca,'LineWidth',1);

set(gcf,'Position',[458,172,639.5,526.5]);

% statistic testing

%% figure 4c d extended data figure 2
% real data
% L = load('data/results/numPattern_results_new/voxel/lh_numPattern_individual_REST1_LR_voxel.mat');
% R = load('data/results/numPattern_results_new/voxel/rh_numPattern_individual_REST1_LR_voxel.mat');

% surrogate data
L = load('data/results/numPattern_results/voxel/surrogate/lh_numPattern_individual_REST1_LR_voxel.mat');
R = load('data/results/numPattern_results/voxel/surrogate/rh_numPattern_individual_REST1_LR_voxel.mat');

patternNames = {'Diverging','Converging','CC-spiral','C-spiral'};


SourceMap = cat(2,L.SourceMap,R.SourceMap);
SinkMap = cat(2,L.SinkMap,R.SinkMap);
CCVortexMap = cat(2,L.CCVortexMap,R.CCVortexMap);
CVortexMap = cat(2,L.CVortexMap,R.CVortexMap);

scale=1;
Size=15;
patternFontsize = 15;
Clim_div=[0 0.08];Clim_curl=[0 0.06];
figure('color','w');
% tiledlayout(2,2,"TileSpacing","compact","Padding","compact");
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);

axes(ha(1));%nexttile;%subplot(2,2,1);
imagesc_brainimg(mean(SourceMap,3),V_mask,0);
% title(patternNames{1},'Fontsize',patternFontsize);
hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim_div);
% set(gca,'outerposition',scale*[0 0.4 0.6 0.6]);

axes(ha(2));%nexttile;%subplot(2,2,2);
imagesc_brainimg(mean(CCVortexMap,3),V_mask,0);
% title(patternNames{3},'Fontsize',patternFontsize);
hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim_curl);
% set(gca,'outerposition',scale*[0.4 0.4 0.6 0.6]);

axes(ha(3));%nexttile;%subplot(2,2,3);
imagesc_brainimg(mean(SinkMap,3),V_mask,0);
% title(patternNames{2},'Fontsize',patternFontsize);
hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim_div);
% set(gca,'outerposition',scale*[0 0 0.6 0.6]);


axes(ha(4));%nexttile;%subplot(2,2,4);
imagesc_brainimg(mean(CVortexMap,3),V_mask,0);
% title(patternNames{4},'Fontsize',patternFontsize);
hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim_curl);
% set(gca,'outerposition',scale*[0.4 0 0.6 0.6]);



colormap(parula)
set(gcf,'position',[36.5,86.5,1342,539]);