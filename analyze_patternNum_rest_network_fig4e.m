clear;clc

addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = 'REST1';
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

%%
L = load(['data/results/numPattern_results_new/network/' hemisphere{1} '_numPattern_individual_' Condition '_LR_network.mat']);
R = load(['data/results/numPattern_results_new/network/' hemisphere{2} '_numPattern_individual_' Condition '_LR_network.mat']);
Source_network = L.Source_network + R.Source_network;
Sink_network = L.Sink_network + R.Sink_network;
CCVortex_network = L.CCVortex_network + R.CCVortex_network;
CVortex_network = L.CVortex_network + R.CVortex_network;
clearvars L R

for i_network = max(label,[],"all"):-1:1
    meanSource(i_network,:) = mean(squeeze(Source_network(i_network,101:1100,:)),1);
end
%%
condition_names = {'Visual', 'Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontalparietal','Default Mode'};

for i_network = max(label,[],"all"):-1:1
    data{i_network} = [reshape(Source_network(i_network,:,:),[],1)...
        reshape(Sink_network(i_network,:,:),[],1)...
        % reshape(CCVortex_network(i_network,:,:),[],1)...
        % reshape(CVortex_network(i_network,:,:),[],1)...
        ];
end

%%
figure;
h = daviolinplot(data,...
    'xtlabels', pattern_names(1:2), ...
    'violin','full',...
    'violinalpha',0.7,...
    'boxcolors','w',...
    'color',parula(7),...
    'outliers',0,...
    'smoothing',1, ...
    'boxspacing',1, ...
    'boxwidth',1);
ylabel('The number of patterns');
% ylim([0 60]);
xl = xlim; xlim([xl(1)-0.05, xl(2)+0.5]);    % make space for the legend
legend([h.ds(1,:)],condition_names);    % add the legend manually
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gcf,'Position',[261,278,739,189.5]);