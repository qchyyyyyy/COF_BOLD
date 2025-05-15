clear;clc

addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = {'WM','RELATIONAL','LANGUAGE','MOTOR'};
PE = 'LR';
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
V_label_L = max(label_L(1:end-1,1:end-1),max(label_L(2:end,1:end-1),max(label_L(1:end-1,2:end),label_L(2:end,2:end))));
V_label_L(V_mask_L==0)=0;
V_label_R = max(label_R(1:end-1,1:end-1),max(label_R(2:end,1:end-1),max(label_R(1:end-1,2:end),label_R(2:end,2:end))));
V_label_R(V_mask_R==0)=0;
V_label = [V_label_L V_label_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R V_label_L V_label_R

pattern_names = {'Diverging', 'Converging', 'CC-Spiral', 'C-Spiral'};
pattern_names_old = {'Source', 'Sink', 'CCVortex', 'CVortex'};

figure;
tiledlayout(length(Condition),1)

for cdn = 1:length(Condition)
    load(['data/results/' Condition{cdn} '_' PE '_SubList.mat']);
    Task.SubList = SubjectList;
    load('data/results/REST1_LR_SubList.mat');
    REST.SubList = SubjectList;
    SubjectList = intersect(Task.SubList,REST.SubList);
    %%
    Subind = find(ismember(Task.SubList,SubjectList));
    Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_' Condition{cdn} '_' PE '.mat']);
    Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_' Condition{cdn} '_' PE '.mat']);
    NumSource = Num_L.NumSource + Num_R.NumSource;
    NumSink = Num_L.NumSink + Num_R.NumSink;
    NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
    NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;

    Task.NumSource = NumSource(Subind,10:end-9);
    Task.NumSink = NumSink(Subind,10:end-9);
    Task.NumCCVortex = NumCCVortex(Subind,10:end-9);
    Task.NumCVortex = NumCVortex(Subind,10:end-9);

    Subind = find(ismember(REST.SubList,SubjectList));
    Num_L = load(['data/results/numPattern_results_new/global/' hemisphere{1} '_numPattern_individual_REST1_LR.mat']);
    Num_R = load(['data/results/numPattern_results_new/global/' hemisphere{2} '_numPattern_individual_REST1_LR.mat']);
    NumSource = Num_L.NumSource + Num_R.NumSource;
    NumSink = Num_L.NumSink + Num_R.NumSink;
    NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
    NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;

    REST.NumSource = NumSource(Subind,101:1100);
    REST.NumSink = NumSink(Subind,101:1100);
    REST.NumCCVortex = NumCCVortex(Subind,101:1100);
    REST.NumCVortex = NumCVortex(Subind,101:1100);


    %%
    condition_names = {Condition{cdn}, 'REST'};

    data{1} = [Task.NumSource(:) Task.NumSink(:) Task.NumCCVortex(:) Task.NumCVortex(:)];
    data{2} = [REST.NumSource(:) REST.NumSink(:) REST.NumCCVortex(:) REST.NumCVortex(:)];

    %%
    nexttile
    h = daviolinplot(data,...
        'xtlabels', pattern_names, ...
        'violin','full',...
        'violinalpha',0.7,...
        'boxcolors','w',...
        'outliers',0,...
        'smoothing',1, ...
        'boxwidth',1.2);
    ylabel('The number of patterns');
    ylim([0 60]);
    xl = xlim; xlim([xl(1)-0.05, xl(2)+0.5]);    % make space for the legend
    legend([h.ds(1,:)],condition_names);    % add the legend manually
    set(gca,'FontSize',8);
    set(gcf,'Position',[261,278,739,189.5]);

    %%

    [h,p(cdn,1),~,stats] = ttest2(REST.NumSource(:),Task.NumSource(:),'Vartype','unequal');t(cdn,1) = stats.tstat;
    [h,p(cdn,2),~,stats] = ttest2(REST.NumSink(:),Task.NumSink(:),'Vartype','unequal');t(cdn,2) = stats.tstat;
    [h,p(cdn,3),~,stats] = ttest2(REST.NumCCVortex(:),Task.NumCCVortex(:),'Vartype','unequal');t(cdn,3) = stats.tstat;
    [h,p(cdn,4),~,stats] = ttest2(REST.NumCVortex(:),Task.NumCVortex(:),'Vartype','unequal')
    t(cdn,4) = stats.tstat

    % Subind2 = 1:391;
    % [h,p(cdn,1),~,stats] = ttest(mean(REST.NumSource(Subind2,:),2),mean(Task.NumSource(Subind2,:),2));t(cdn,1) = stats.tstat;
    % [h,p(cdn,2),~,stats] = ttest(mean(REST.NumSink(Subind2,:),2),mean(Task.NumSink(Subind2,:),2));t(cdn,2) = stats.tstat;
    % [h,p(cdn,3),~,stats] = ttest(mean(REST.NumCCVortex(Subind2,:),2),mean(Task.NumCCVortex(Subind2,:),2));t(cdn,3) = stats.tstat;
    % [h,p(cdn,4),~,stats] = ttest(mean(REST.NumCVortex(Subind2,:),2),mean(Task.NumCVortex(Subind2,:),2))
    % t(cdn,4) = stats.tstat

end
set(gcf,'Position',[261,0,739,189.5*length(Condition)]);