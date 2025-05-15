clear;clc
addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = 'WM';% 'WM' 'RELATIONAL'
switch Condition
    case 'WM'
        TaskName = 'Working Memory';
    case 'RELATIONAL'
        TaskName = 'Relational Processing';
end

PE = 'LR';% 'RL'
switch PE
    case 'RL'
        Run = '1';
    case 'LR'
        Run = '2';
end

[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
load(['data/results/info_task/' Condition '_' PE '_info_new.mat']);
load(['data/results/' Condition '_' PE '_Sublist_new.mat']);

% global cortex level
if length(hemisphere)==2
    Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition '_' PE '.mat']);
    Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_' Condition '_' PE '.mat']);
    NumSource = Num_L.NumSource + Num_R.NumSource;
    NumSink = Num_L.NumSink + Num_R.NumSink;
    NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
    NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;
else
    load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition '_' PE '.mat']);
end
%%
if strcmp(Condition,'WM')
    Acc(Acc==0)=nan;
end
Subind=find(~isnan(Acc(:,6)));% 1:size(NumSource,1);
T = 5:size(NumSource,2)-4;
T(T>size(NumSource,2))=[];
NumSource = NumSource(Subind,T);
NumSink = NumSink(Subind,T);
NumCCVortex = NumCCVortex(Subind,T);
NumCVortex = NumCVortex(Subind,T);

Acc = Acc(Subind,~isnan(Acc(1,:)));
[transformed_data, lambda] = boxcox(mean(Acc,2,'omitnan'));
% Acc = log(Acc ./ (1 - Acc));
RT = RT(Subind,~isnan(RT(1,:)));
Threshold = 0;

SubjectList = SubjectList(Subind);
%%
behav_info = readtable('data/fMRI/unrestricted_qchyyyyyy_11_22_2021_0_28_6.csv');
behav_info = behav_info(ismember(behav_info.Subject,str2num(cell2mat(SubjectList))),:);

sub_gender = behav_info.Gender;
sub_gender_bi = zeros(size(sub_gender));
sub_gender_bi(cellfun(@(x) x=='M',sub_gender))=1;
sub_age = behav_info.Age;
sub_age = cellfun(@(x) split(x,'-'), sub_age,'UniformOutput',false);
[sub_age{cell2mat(cellfun(@(x) length(x)<2,sub_age,'UniformOutput',false))}] = deal({'36';'36'});
sub_age = cell2mat(cellfun(@(x) (str2double(x{1})+str2double(x{2}))/2, sub_age,'UniformOutput',false));
%%
pattern_names_old = {'Source','Sink','CCVortex','CVortex'};
pattern_names = {'diverging','converging','CC-spiral','C-spiral'};

figure;
tt = tiledlayout(2,5,"TileSpacing","compact","Padding",'compact');
xlabel(tt,'Number of Patterns','FontSize',15);

title(tt,[TaskName ' Run ' Run],'FontSize',18);

%
YYLim = [mean(transformed_data)-3*std(transformed_data) mean(transformed_data)+3*std(transformed_data)];

for i_pattern = 1:4
    nexttile
    eval(['x=mean(Num' pattern_names_old{i_pattern} ',2);'])
    Y=transformed_data;
    X_design = [x sub_age sub_gender_bi];

    [b,stats] = robustfit(X_design,Y);
    t(i_pattern,1) = stats.t(2);
    P(i_pattern,1) = stats.p(2);

    XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
    x_fit = linspace(XXLim(1), XXLim(2), 100)';
    X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
    y_fit = [ones(size(x_fit)), X_pred] * b;

    h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.5);alpha(h3,0.5);
    hold on
    h2=plot(x_fit,y_fit,'b','linewidth',2);
    xlim(XXLim);ylim(YYLim);
    title(pattern_names{i_pattern},'FontSize',15);
    if i_pattern==1;ylabel('Accuracy (boxcox)','FontSize',15);end

    text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2))]},'FontSize',10);
    set(gca,'LineWidth',1);
end

nexttile
x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);
Y=transformed_data;
X_design = [x sub_age sub_gender_bi];

[b,stats] = robustfit(X_design,Y);
t(5,1) = stats.t(2);
P(5,1) = stats.p(2);

XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
x_fit = linspace(XXLim(1), XXLim(2), 100)';
X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
y_fit = [ones(size(x_fit)), X_pred] * b;

h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.5);alpha(h3,0.5);
hold on
h2=plot(x_fit,y_fit,'b','linewidth',2);
xlim(XXLim);ylim(YYLim);
title('all','FontSize',15);
text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2))]},'FontSize',10);
set(gca,'LineWidth',1);
%%

for i_pattern = 1:4
    nexttile
    eval(['x=mean(Num' pattern_names_old{i_pattern} ',2);'])
    Y=mean(RT,2,'omitnan');
    X_design = [x sub_age sub_gender_bi];

    [b,stats] = robustfit(X_design,Y);
    t(i_pattern,2) = stats.t(2);
    P(i_pattern,2) = stats.p(2);

    XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
    YYLim = [mean(Y)-3*std(Y) mean(Y)+3*std(Y)];

    x_fit = linspace(XXLim(1), XXLim(2), 100)';
    X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
    y_fit = [ones(size(x_fit)), X_pred] * b;

    h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.5);alpha(h3,0.5);
    hold on
    h2=plot(x_fit,y_fit,'b','linewidth',2);
    xlim(XXLim);ylim(YYLim);
    title(pattern_names{i_pattern},'FontSize',15);
    if i_pattern==1;ylabel('Reaction time','FontSize',15);end

    % ylabel('Surrogate','FontSize',10);
    text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2))]},'FontSize',10);
    set(gca,'LineWidth',1);
end

nexttile
x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);
Y=mean(RT,2,'omitnan');
X_design = [x sub_age sub_gender_bi];

[b,stats] = robustfit(X_design,Y);
t(5,2) = stats.t(2);
P(5,2) = stats.p(2);

XXLim = [mean(x)-3*std(x), mean(x)+3*std(x)];
YYLim = [mean(Y)-3*std(Y) mean(Y)+3*std(Y)];

x_fit = linspace(XXLim(1), XXLim(2), 100)';
X_pred = [x_fit, mean(sub_age)*ones(size(x_fit)), mean(sub_gender_bi)*ones(size(x_fit))];
y_fit = [ones(size(x_fit)), X_pred] * b;

h3=scatter(x,Y,[],[0.6 0.6 0.6],"o",'filled','LineWidth',0.5);alpha(h3,0.5);
hold on
h2=plot(x_fit,y_fit,'b','linewidth',2);
xlim(XXLim);ylim(YYLim);
title('all','FontSize',15);
text(0.95*XXLim(1)+0.05*XXLim(2),0.95*YYLim(2)+0.05*YYLim(1),{['P = ' num2str(stats.p(2))]},'FontSize',10);
set(gca,'LineWidth',1);
set(gcf,'Position',[2.5,27,1408,331*2]);