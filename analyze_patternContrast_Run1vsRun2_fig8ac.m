clear;clc

addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = 'WM';
PE = {'RL','LR'};
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
%%

load(['data/results/info_task/' Condition '_' PE{1} '_info.mat']);
Subind = [1:135 137:216];
RL.Acc = mean(Acc(Subind,~isnan(Acc(1,:))),2);
RL.RT = mean(RT(Subind,~isnan(Acc(1,:))),2);
load(['data/results/info_task/' Condition '_' PE{2} '_info.mat']);
Subind = 1:215;
LR.Acc = mean(Acc(Subind,~isnan(Acc(1,:))),2);
LR.RT = mean(RT(Subind,~isnan(Acc(1,:))),2);

%% figure 8a

% Accuracy
[h,p,~,stat] = ttest(RL.Acc(~isnan(RL.Acc) & ~isnan(LR.Acc)),LR.Acc(~isnan(RL.Acc) & ~isnan(LR.Acc)))

g1 = repmat({'Run1'},length(find(~isnan(RL.Acc) & ~isnan(LR.Acc))),1);
g2 = repmat({'Run2'},length(find(~isnan(RL.Acc) & ~isnan(LR.Acc))),1);
figure;
set(gca,'Linewidth',1.1);
bb=boxplot([RL.Acc(~isnan(RL.Acc) & ~isnan(LR.Acc));LR.Acc(~isnan(RL.Acc) & ~isnan(LR.Acc))],[g1;g2],"OutlierSize",0.01);
set(bb,'LineWidth',1);
hold on
sigline([1,2],1,p,0.1);
ylabel('Accuracy');
ylim([0.6 1.2]);
set(gca,'Linewidth',1);
set(gcf,'Position',[896,255,167.5,230])

% Reaction time
[h,p,~,stat] = ttest(RL.RT(~isnan(RL.RT) & ~isnan(LR.RT)),LR.Acc(~isnan(RL.RT) & ~isnan(LR.RT)))

g1 = repmat({'Run1'},length(find(~isnan(RL.RT) & ~isnan(LR.RT))),1);
g2 = repmat({'Run2'},length(find(~isnan(RL.RT) & ~isnan(LR.RT))),1);
figure;
bb = boxplot([RL.RT(~isnan(RL.RT) & ~isnan(LR.RT));LR.RT(~isnan(RL.RT) & ~isnan(LR.RT))],[g1;g2],"OutlierSize",0.01);
set(bb,'LineWidth',1);
hold on
sigline([1,2],1300,p,0.1);
ylabel('Reaction time');
ylim([500 1600]);
set(gca,'Linewidth',1);
set(gcf,'Position',[896,255,167.5,230])
%% voxel-wise PQ Map
Subind = [1:135 137:216];

L = load(['data/results/numPattern_results/voxel/' hemisphere{1} '_numPattern_individual_' Condition '_' PE{1} '_voxel.mat']);
R = load(['data/results/numPattern_results/voxel/' hemisphere{2} '_numPattern_individual_' Condition '_' PE{1} '_voxel.mat']);
RL.SourceMap = cat(2,L.SourceMap(:,:,Subind),R.SourceMap(:,:,Subind));
RL.SinkMap = cat(2,L.SinkMap(:,:,Subind),R.SinkMap(:,:,Subind));
RL.CCVortexMap = cat(2,L.CCVortexMap(:,:,Subind),R.CCVortexMap(:,:,Subind));
RL.CVortexMap = cat(2,L.CVortexMap(:,:,Subind),R.CVortexMap(:,:,Subind));

Subind = 1:215;

L = load(['data/results/numPattern_results/voxel/' hemisphere{1} '_numPattern_individual_' Condition '_' PE{2} '_voxel.mat']);
R = load(['data/results/numPattern_results/voxel/' hemisphere{2} '_numPattern_individual_' Condition '_' PE{2} '_voxel.mat']);
LR.SourceMap = cat(2,L.SourceMap(:,:,Subind),R.SourceMap(:,:,Subind));
LR.SinkMap = cat(2,L.SinkMap(:,:,Subind),R.SinkMap(:,:,Subind));
LR.CCVortexMap = cat(2,L.CCVortexMap(:,:,Subind),R.CCVortexMap(:,:,Subind));
LR.CVortexMap = cat(2,L.CVortexMap(:,:,Subind),R.CVortexMap(:,:,Subind));

%% figure 8c mass univariate analysis (paired t-test, run 1 vs run 2)

pattern_names = {'diverging','converging','CC-spiral','C-spiral'};

scale=1;
Size=5;
patternFontsize = 16;
figure('color','w');
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);
Clim = [-10 10];
axes(ha(1));
[~,P,~,stats] = ttest(map2vec(LR.SourceMap,V_mask),map2vec(RL.SourceMap,V_mask),"Dim",2);
    Tmap = vec2map(stats.tstat,V_mask);Tmap(V_mask==0)=0;
    P = FDR_correction(P);
    Pmap = vec2map(P,V_mask);Pmap(V_mask==0)=nan;
% imagesc_brainimg(-log10(Pmap),Tmap>0,0);
imagesc_brainimg(Tmap,Pmap<0.05,0);
title(pattern_names{1},'Fontsize',patternFontsize);hold on;scatter_boundary(oi_x,oi_y,Size);
clim(Clim);
% set(gca,'outerposition',scale*[0 0 0.25 1]);

axes(ha(3));
[~,P,~,stats] = ttest(map2vec(LR.SinkMap,V_mask),map2vec(RL.SinkMap,V_mask),"Dim",2);
    Tmap = vec2map(stats.tstat,V_mask);Tmap(V_mask==0)=0;
        P = FDR_correction(P);
    Pmap = vec2map(P,V_mask);Pmap(V_mask==0)=nan;
% imagesc_brainimg(-log10(Pmap),Tmap>0,0);
imagesc_brainimg(Tmap,Pmap<0.05,0);
title(pattern_names{2},'Fontsize',patternFontsize);hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim);
% set(gca,'outerposition',scale*[0.25 0 0.25 1]);

axes(ha(2));
[~,P,~,stats] = ttest(map2vec(LR.CCVortexMap,V_mask),map2vec(RL.CCVortexMap,V_mask),"Dim",2);
    Tmap = vec2map(stats.tstat,V_mask);Tmap(V_mask==0)=0;
        P = FDR_correction(P);
    Pmap = vec2map(P,V_mask);Pmap(V_mask==0)=nan;
% imagesc_brainimg(-log10(Pmap),Tmap>0,0);
imagesc_brainimg(Tmap,Pmap<0.05,0);
title(pattern_names{3},'Fontsize',patternFontsize);hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim);
% set(gca,'outerposition',scale*[0.5 0 0.25 1]);

axes(ha(4));
[~,P,~,stats] = ttest(map2vec(LR.CVortexMap,V_mask),map2vec(RL.CVortexMap,V_mask),"Dim",2);
    Tmap = vec2map(stats.tstat,V_mask);Tmap(V_mask==0)=0;
        P = FDR_correction(P);
    Pmap = vec2map(P,V_mask);Pmap(V_mask==0)=nan;
% imagesc_brainimg(-log10(Pmap),Tmap>0,0);
imagesc_brainimg(Tmap,Pmap<0.05,0);
title(pattern_names{4},'Fontsize',patternFontsize);hold on;scatter_boundary(oi_x,oi_y,Size);clim(Clim);
% set(gca,'outerposition',scale*[0.75 0 0.25 1]);

% colormap(gcf,hot);
colormap(gcf,bluewhitered);
set(gcf,'position',[365.5,144.5,1062.5,480]);
%% pair test t map 3d
HCPTemplatePath = 'data/templates/HCP_S1200_Atlas_Z4_pkXDZ/';
patternFontsize = 12;
Clim = [-10 10];
pattern_names_old_resort = {'Source','CCVortex','Sink','CVortex'};
pattern_names_resort = {'diverging','CC-spiral','converging','C-spiral'};

figure('color','w');
ha = tight_subplot(4,4,[.01 .01],[.01 .01],[.01 .01]);
for pa = 1:4
    % suptitle(ha(),patternName{pa},'Fontsize',patternFontsize);
    
    eval(['[~,P,~,stats] = ttest(map2vec(LR.' pattern_names_old_resort{pa} 'Map,V_mask),map2vec(RL.' pattern_names_old_resort{pa} 'Map,V_mask),"Dim",2);']);
    Tmap = vec2map(stats.tstat,V_mask);Tmap(V_mask==0)=0;
    P = FDR_correction(P);
    Pmap = vec2map(P,V_mask);Pmap(V_mask==0)=nan;
    Tmap(Pmap>=0.05)=0;
    axes(ha(2*pa-1+4*(pa>2)));
    patch_3dbrainimg(Tmap(:,1:250),...
        [HCPTemplatePath 'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'],...
        [HCPTemplatePath 'S1200.L.flat.32k_fs_LR.surf.gii'],Clim);
    camlight('headlight')
    material dull
    axis off
    axis image

    text(-66,-62,86,pattern_names_resort{pa},"FontSize",16);

    axes(ha(2*pa+4*(pa>2)));
    patch_3dbrainimg(Tmap(:,251:end),...
        [HCPTemplatePath 'S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'],...
        [HCPTemplatePath 'S1200.R.flat.32k_fs_LR.surf.gii'],Clim);[caz,cel] = view();view(-caz,cel);
        camlight('headlight')
    material dull
    axis off
    axis image
    axes(ha(4+2*pa-1+4*(pa>2)));
    patch_3dbrainimg(Tmap(:,1:250),...
        [HCPTemplatePath 'S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'],...
        [HCPTemplatePath 'S1200.L.flat.32k_fs_LR.surf.gii'],Clim);[caz,cel] = view();view(-caz,cel);
        camlight('headlight')
    material dull
    axis off
    axis image
    axes(ha(4++2*pa+4*(pa>2)));
    patch_3dbrainimg(Tmap(:,251:end),...
        [HCPTemplatePath 'S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'],...
        [HCPTemplatePath 'S1200.R.flat.32k_fs_LR.surf.gii'],Clim);
    camlight('headlight')
    material dull
    axis off
    axis image
end

colormap(gcf,[[0.5 0.5 0.5];bluewhitered]);
set(gcf,'position',[297.5,80.5,896,702]);