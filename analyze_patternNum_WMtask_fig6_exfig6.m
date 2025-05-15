clear;clc

addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = 'WM';
PE = 'LR';
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

pattern_names = {'diverging','converging','CC-spiral','C-spiral'};
%% figure 5a
load(['data/results/info_task/' Condition '_' PE '_info.mat']);
Block_name = {'0_Back_Body','0_Back_Face',...
    '0_Back_Place','0_Back_Tools',...
    '2_Back_Body','2_Back_Face',...
    '2_Back_Place','2_Back_Tools'};

% Subind = 1:215;

Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition '_' PE '.mat']);
Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_' Condition '_' PE '.mat']);
NumSource = Num_L.NumSource + Num_R.NumSource;
NumSink = Num_L.NumSink + Num_R.NumSink;
NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;

NumSource(BadInd,:)=[];
NumSink(BadInd,:)=[];
NumCCVortex(BadInd,:)=[];
NumCVortex(BadInd,:)=[];

%%
load(['data/results/' Condition '_' PE '_SubList.mat']);
Task.SubjectList = SubjectList;
load(['data/results/REST1_' PE '_SubList.mat']);
REST.SubjectList = SubjectList;
SubjectList = intersect(Task.SubjectList,REST.SubjectList);
subidx_task = find(ismember(Task.SubjectList,SubjectList));
subidx_rest = find(ismember(REST.SubjectList,SubjectList));

Task.NumSource = NumSource(subidx_task,10:end-9);
Task.NumSink = NumSink(subidx_task,10:end-9);
Task.NumCCVortex = NumCCVortex(subidx_task,10:end-9);
Task.NumCVortex = NumCVortex(subidx_task,10:end-9);

Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_REST1_LR.mat']);
Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_REST1_LR.mat']);
REST.NumSource = Num_L.NumSource + Num_R.NumSource;
REST.NumSink = Num_L.NumSink + Num_R.NumSink;
REST.NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
REST.NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;clearvars Num_L Num_R


REST.NumSource = REST.NumSource(subidx_rest,101:1100);
REST.NumSink = REST.NumSink(subidx_rest,101:1100);
REST.NumCCVortex = REST.NumCCVortex(subidx_rest,101:1100);
REST.NumCVortex = REST.NumCVortex(subidx_rest,101:1100);

%%
condition_names = {Condition, 'REST'};
% pattern_names = {'Source', 'Sink', 'CC-spiral', 'C-spiral'};

data{1} = [Task.NumSource(:) Task.NumSink(:) Task.NumCCVortex(:) Task.NumCVortex(:)];
data{2} = [REST.NumSource(:) REST.NumSink(:) REST.NumCCVortex(:) REST.NumCVortex(:)];

[~,p(1),~,stats] = ttest2(Task.NumSource(:),REST.NumSource(:),'Vartype','unequal');t(1) = stats.tstat;
[~,p(2),~,stats] = ttest2(Task.NumSink(:),REST.NumSink(:),'Vartype','unequal');t(2) = stats.tstat;
[~,p(3),~,stats] = ttest2(Task.NumCCVortex(:),REST.NumCCVortex(:),'Vartype','unequal');t(3) = stats.tstat;
[~,p(4),~,stats] = ttest2(Task.NumCVortex(:),REST.NumCVortex(:),'Vartype','unequal')
t(4) = stats.tstat

% [~,p(1),~,stats] = ttest(mean(Task.NumSource,2),mean(REST.NumSource,2));t(1) = stats.tstat;
% [~,p(2),~,stats] = ttest(mean(Task.NumSink,2),mean(REST.NumSink,2));t(2) = stats.tstat;
% [~,p(3),~,stats] = ttest(mean(Task.NumCCVortex,2),mean(REST.NumCCVortex,2));t(3) = stats.tstat;
% [~,p(4),~,stats] = ttest(mean(Task.NumCVortex,2),mean(REST.NumCVortex,2))
% t(4) = stats.tstat
%%
figure;
subplot(7,1,[1 2]);
h = daviolinplot(data,...
    'xtlabels', pattern_names, ...
    'violin','full',...
    'violinalpha',0.7,...
    'boxcolors','w',...
    'outliers',0,...
    'smoothing',1, ...
    'boxwidth',1.2); 
% ylabel('The number of patterns');
ylim([0 80]);
xl = xlim; xlim([xl(1)-0.05, xl(2)+0.5]);    % make space for the legend
ll = legend([h.ds(1,:)],condition_names);    % add the legend manually
ll.AutoUpdate='off';
Int = 0.2;
Idx = [1 1];
% hold on
% sigline(1,70,p(1),Int,Idx);
% hold on
% sigline(2,70,p(2),Int,Idx);
% hold on
% sigline(3,70,p(3),Int,Idx);
% hold on
% sigline(4,70,p(4),Int,Idx);


set(gca,'FontSize',8);
ax = gca;
ax.XAxis.FontSize = 14;
set(gca,'LineWidth',1);

%% figure 5c
Delay = 18;

Events = zeros(404,1);

for bb = 1:8
    idx = ceil(meanStimType(bb,1)-0.5):fix(meanStimType(bb,2)-0.5);
    idx = idx+Delay;
    idx(idx>404)=[];idx(idx<1)=[];
    Events(idx) = meanStimType(bb,3);
end

WM_0b.NumSource = NumSource(:,Events<5 & Events>0);
WM_0b.NumSink = NumSink(:,Events<5 & Events>0);
WM_0b.NumCCVortex = NumCCVortex(:,Events<5 & Events>0);
WM_0b.NumCVortex = NumCVortex(:,Events<5 & Events>0);

WM_2b.NumSource = NumSource(:,Events>4);
WM_2b.NumSink = NumSink(:,Events>4);
WM_2b.NumCCVortex = NumCCVortex(:,Events>4);
WM_2b.NumCVortex = NumCVortex(:,Events>4);

Fixation.NumSource = NumSource(:,Events==0);
Fixation.NumSink = NumSink(:,Events==0);
Fixation.NumCCVortex = NumCCVortex(:,Events==0);
Fixation.NumCVortex = NumCVortex(:,Events==0);


condition_names = {'WM_2b', 'WM_0b', 'Fixation'};
% pattern_names = {'Source', 'Sink', 'CC-spiral', 'C-spiral'};

data{1} = [WM_2b.NumSource(:) WM_2b.NumSink(:) WM_2b.NumCCVortex(:) WM_2b.NumCVortex(:)];
data{2} = [WM_0b.NumSource(:) WM_0b.NumSink(:) WM_0b.NumCCVortex(:) WM_0b.NumCVortex(:)];
data{3} = [Fixation.NumSource(:) Fixation.NumSink(:) Fixation.NumCCVortex(:) Fixation.NumCVortex(:)];

subplot(7,1,[6 7]);

h = daviolinplot(data,...
    'xtlabels', pattern_names, ...
    'violin','full',...
    'violinalpha',0.7,...
    'boxcolors','w',...
    'outliers',0,...
    'smoothing',1, ...
    'boxwidth',1.2); 

ylim([0 80]);
xl = xlim; xlim([xl(1)-0.05, xl(2)+0.5]);    % make space for the legend
ll = legend([h.ds(1,:)],condition_names,'Interpreter','none');    % add the legend manually
ll.AutoUpdate='off';

% Int = 0.2;
% Idx = [1 0 1];
% hold on
% sigline(1,80,p(1),Int,Idx);
% hold on
% sigline(2,80,p(2),Int,Idx);
% hold on
% sigline(3,80,p(3),Int,Idx);
% hold on
% sigline(4,80,p(4),Int,Idx);
% hold on
% Int = 0.2;
% Idx = [0 1 1];
% hold on
% sigline(1,70,p(1),Int,Idx);
% hold on
% sigline(2,70,p(2),Int,Idx);
% hold on
% sigline(3,70,p(3),Int,Idx);
% hold on
% sigline(4,70,p(4),Int,Idx);
% hold on
% Int = 0.2;
% Idx = [1 1 0];
% hold on
% sigline(1,60,p(1),Int,Idx);
% hold on
% sigline(2,60,p(2),Int,Idx);
% hold on
% sigline(3,60,p(3),Int,Idx);
% hold on
% sigline(4,60,p(4),Int,Idx);

set(gca,'FontSize',8);
ax = gca;
ax.XAxis.FontSize = 14;
set(gca,'LineWidth',1);
%% figure 5b (temporal fluctuation of pattern quantity)
Time = 5:400;
Time(Time>404)=[];
Time(Time<1)=[];
NumSource = NumSource(:,Time);
NumSink = NumSink(:,Time);
NumCCVortex = NumCCVortex(:,Time);
NumCVortex = NumCVortex(:,Time);

YY = [0 60];
CC = [0.92 0.69 0.12;...
    0.85 0.32 0.10;...
    0.5 0.8 0.3;
    0 0.44 0.74];
Delay = 18;

subplot(7,1,[3 4 5]);
set(gca,'FontSize',8);

for i = 1:length(meanStimType)
    XX = [meanStimType(i,1)+Delay meanStimType(i,2)+Delay meanStimType(i,2)+Delay meanStimType(i,1)+Delay];
    pp = patch(XX,[YY(1)+0.5 YY(1)+0.5 YY(2)-0.5 YY(2)-0.5],[0.96 0.96 0.96],'edgecolor','none');hold on
end

p1 = plot(Time,mean(NumSource,1),'Color',CC(1,:),'LineWidth',2);hold on
p2 = plot(Time,mean(NumSink,1),'Color',CC(2,:),'LineWidth',2);hold on
p3 = plot(Time,mean(NumCVortex,1),'Color',CC(3,:),'LineWidth',2);hold on
p4 = plot(Time,mean(NumCCVortex,1),'Color',CC(4,:),'LineWidth',2);hold on

ll = legend([pp,p1,p2,p3,p4],'Task block',pattern_names{1},pattern_names{2},pattern_names{3},pattern_names{4});
set(ll,'AutoUpdate','off');

pp = patch([Time flip(Time)],[mean(NumSource,1)+std(NumSource,[],1) flip(mean(NumSource,1)-std(NumSource,[],1))],...
    CC(1,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumSink,1)+std(NumSink,[],1) flip(mean(NumSink,1)-std(NumSink,[],1))],...
    CC(2,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCCVortex,1)+std(NumCCVortex,[],1) flip(mean(NumCCVortex,1)-std(NumCCVortex,[],1))],...
    CC(3,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCVortex,1)+std(NumCVortex,[],1) flip(mean(NumCVortex,1)-std(NumCVortex,[],1))],...
    CC(4,:),'edgecolor','none');hold on
alpha(pp,0.2);



xlim([0 404]);
ylim(YY);

xlabel('Time (TR)','FontSize',14);
yy = ylabel(gca,'The number of patterns on whole cortex per frame','FontSize',16);
set(yy,'position',[-29.304508499630458,28.60352422907491,1.4e-14])
set(gca,'LineWidth',1);

set(gcf,"Position",[341.5,-196.5,873,904]);
%%
Delay = 0;

Events = zeros(404,1);

for bb = 1:8
    idx = ceil(meanStimType(bb,1)-0.5):fix(meanStimType(bb,2)-0.5);
    idx = idx+Delay;
    idx(idx>404)=[];idx(idx<1)=[];
    Events(idx) = meanStimType(bb,3);
end
Events = Events(Time);
Events(Events<=4 & Events>=1)=1;
Events(Events>=5)=2;

%%
NumPattern = NumSource;

[R,lag] = crosscorr(mean(NumPattern,1),Events);
R(3)
figure;crosscorr(mean(NumPattern,1),Events,NumLags=30);

for sub = size(NumPattern):-1:1
    [R,lag] = crosscorr(NumPattern(sub,:)',Events);
    R_sub(sub) = R(3);
end
mean(R_sub)
std(R_sub)
%% extended data figure 6
figure;
subplot(2,1,1);
Delay = 0;title('Original','FontSize',14)
set(gca,'FontSize',8);

for i = 1:length(meanStimType)
    XX = [meanStimType(i,1)+Delay meanStimType(i,2)+Delay meanStimType(i,2)+Delay meanStimType(i,1)+Delay];
    pp = patch(XX,[YY(1)+0.5 YY(1)+0.5 YY(2)-0.5 YY(2)-0.5],[0.96 0.96 0.96],'edgecolor','none');hold on
end

p1 = plot(Time,mean(NumSource,1),'Color',CC(1,:),'LineWidth',2);hold on
p2 = plot(Time,mean(NumSink,1),'Color',CC(2,:),'LineWidth',2);hold on
p3 = plot(Time,mean(NumCVortex,1),'Color',CC(3,:),'LineWidth',2);hold on
p4 = plot(Time,mean(NumCCVortex,1),'Color',CC(4,:),'LineWidth',2);hold on

ll = legend([pp,p1,p2,p3,p4],'Task block',pattern_names{1},pattern_names{2},pattern_names{3},pattern_names{4});
set(ll,'AutoUpdate','off');

pp = patch([Time flip(Time)],[mean(NumSource,1)+std(NumSource,[],1) flip(mean(NumSource,1)-std(NumSource,[],1))],...
    CC(1,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumSink,1)+std(NumSink,[],1) flip(mean(NumSink,1)-std(NumSink,[],1))],...
    CC(2,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCCVortex,1)+std(NumCCVortex,[],1) flip(mean(NumCCVortex,1)-std(NumCCVortex,[],1))],...
    CC(3,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCVortex,1)+std(NumCVortex,[],1) flip(mean(NumCVortex,1)-std(NumCVortex,[],1))],...
    CC(4,:),'edgecolor','none');hold on
alpha(pp,0.2);

xlim([0 404]);
ylim(YY);

subplot(2,1,2);
Delay = 18;title('Time shift of 18 TRs','FontSize',14)
set(gca,'FontSize',8);

for i = 1:length(meanStimType)
    XX = [meanStimType(i,1)+Delay meanStimType(i,2)+Delay meanStimType(i,2)+Delay meanStimType(i,1)+Delay];
    pp = patch(XX,[YY(1)+0.5 YY(1)+0.5 YY(2)-0.5 YY(2)-0.5],[0.96 0.96 0.96],'edgecolor','none');hold on
end

p1 = plot(Time,mean(NumSource,1),'Color',CC(1,:),'LineWidth',2);hold on
p2 = plot(Time,mean(NumSink,1),'Color',CC(2,:),'LineWidth',2);hold on
p3 = plot(Time,mean(NumCVortex,1),'Color',CC(3,:),'LineWidth',2);hold on
p4 = plot(Time,mean(NumCCVortex,1),'Color',CC(4,:),'LineWidth',2);hold on

ll = legend([pp,p1,p2,p3,p4],'Task block',pattern_names{1},pattern_names{2},pattern_names{3},pattern_names{4});
set(ll,'AutoUpdate','off');

pp = patch([Time flip(Time)],[mean(NumSource,1)+std(NumSource,[],1) flip(mean(NumSource,1)-std(NumSource,[],1))],...
    CC(1,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumSink,1)+std(NumSink,[],1) flip(mean(NumSink,1)-std(NumSink,[],1))],...
    CC(2,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCCVortex,1)+std(NumCCVortex,[],1) flip(mean(NumCCVortex,1)-std(NumCCVortex,[],1))],...
    CC(3,:),'edgecolor','none');hold on
alpha(pp,0.2);
pp = patch([Time flip(Time)],[mean(NumCVortex,1)+std(NumCVortex,[],1) flip(mean(NumCVortex,1)-std(NumCVortex,[],1))],...
    CC(4,:),'edgecolor','none');hold on
alpha(pp,0.2);

xlim([0 404]);
ylim(YY);

xlabel('Time (TR)','FontSize',14);
yy = ylabel(gca,'The number of patterns on whole cortex per frame','FontSize',14);
set(yy,'position',[-28.701404286770128,70.59319286871961,1.4e-14])
set(gca,'LineWidth',1);
set(gcf,"Position",[137.5,82.5,918.5,711]);