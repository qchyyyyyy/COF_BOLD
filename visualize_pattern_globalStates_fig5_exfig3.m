clear;
addpath(genpath('functions'));

hemisphere = {'lh','rh'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks'
[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere{1});
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere{2});
PE = 'LR';
mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

OutputPath = 'data/results/';

%% Load data
load([OutputPath 'clustering_results/LR_patterns_Kmeans_correlation_K=20.mat']);

%% figure 5

figure('color','w');
ha = tight_subplot(2,2,[.01 .01],[.01 .01],[.01 .01]);
selected_mode = [5 9 11 12];
pattern_order = [1 3 2 4];
for mm = 1:length(selected_mode)
    figure('color','w');
    set(gcf, 'Renderer', 'painters');
    set(gcf,'position',[102.5,71.5,1169/2,632/2]);

    DivMap = vec2map(C(selected_mode(mm),1:length(find(V_mask(:)==1)))',V_mask);
    CurlMap = vec2map(C(selected_mode(mm),length(find(V_mask(:)==1))+1:end)',V_mask);

    DivMap(V_mask==0)=nan;DivMap = (DivMap-mean(DivMap,'all','omitnan'))./std(DivMap,[],'all','omitnan');
    CurlMap(V_mask==0)=nan;CurlMap = (CurlMap-mean(CurlMap,'all','omitnan'))./std(CurlMap,[],'all','omitnan');
    
    [Ux,Uy] = recon_velfield(DivMap,CurlMap,V_mask);
    
    % axes(ha(mm));
    ax3 = gca;
    imagesc_brainimg(mask,1-mask,0,[0 1]);colormap(ax3,[[1 1 1]; jet]);
    ax3.Color='k';

    ax1=axes;
    imagesc_brainimg(CurlMap,V_mask,0);
    alpha(abs(CurlMap)/3);
    colormap(ax1,winter);
    ax1.Color='none';ax1.XColor='none';ax1.YColor='none';
    axis([0 250 0 200]);axis equal;axis tight;
    
    ax2=axes;
    imagesc_brainimg(DivMap,V_mask,0);
    alpha(abs(DivMap)/3);
    colormap(ax2,autumn);
    ax2.Color='none';ax2.XColor='none';ax2.YColor='none';
    axis([0 250 0 200]);axis equal;axis tight;
    
    hold on;
    scatter_boundary(oi_x,oi_y,20);

    hold on
    q = streamslice(Ux,Uy,4);
    set(q,'Color',[0.9 0.9 0.9]);
    hold off
    % tt = title(['motif ' num2str(pattern_order(mm)) ],'FontSize',15);
    set(ax1,'position',ax3.Position);
    set(ax2,'position',ax3.Position);
end

%% extended data figure 3
figure('color','w');
set(gcf,'Units','centimeters');
set(gcf,'Position',[-1,-1,18,22.5]);
% [~,idx_states] = sort(mean_state,'descend');
idx_states = 1:20;
X = 2;Y = 10;Px=0.7;Py=0.7;

ha = tight_subplot(7,3,[.01 .01],[.01 .01],[.01 .01]);
for m = 1:20
    mm = idx_states(m);
    DivMap = vec2map(C(mm,1:length(find(V_mask(:)==1)))',V_mask);
    CurlMap = vec2map(C(mm,length(find(V_mask(:)==1))+1:end)',V_mask);

    DivMap(V_mask==0)=nan;DivMap = (DivMap-mean(DivMap,'all','omitnan'))./std(DivMap,[],'all','omitnan');
    CurlMap(V_mask==0)=nan;CurlMap = (CurlMap-mean(CurlMap,'all','omitnan'))./std(CurlMap,[],'all','omitnan');
    
    axes(ha(mm));
    ax3 = gca;
    imagesc_brainimg(V_mask,1-V_mask,0,[0 1]);colormap(ax3,[[1 1 1]; jet]);
    ax3.Color='k';
    set(gca,'FontName','Arial','FontUnits','points','FontSize',6);

    ax1=axes;
    imagesc_brainimg(CurlMap,V_mask,0);
    alpha(abs(CurlMap)/3);
    colormap(ax1,winter);
    ax1.Color='none';ax1.XColor='none';ax1.YColor='none';
    axis([0 250 0 200]);axis equal;axis tight;
    set(ax1,'position',ax3.Position);
    set(gca,'FontName','Arial','FontUnits','points','FontSize',6);

    ax2=axes;
    imagesc_brainimg(DivMap,V_mask,0);
    alpha(abs(DivMap)/3);
    colormap(ax2,autumn);
    ax2.Color='none';ax2.XColor='none';ax2.YColor='none';
    axis([0 250 0 200]);axis equal;axis tight;
    set(ax2,'position',ax3.Position);

    hold on;
    ss = scatter(oi_x,oi_y,0.01,'.');
    alpha(ss,0.1);
    ss.MarkerEdgeColor = [0.95 0.95 0.95];
    tt = title(['Motif ' num2str(m)]);
    % tt = title(['motif ' num2str(m) ' (' num2str(perc_states_group(mm)*100) '%)']);
    tt.Color='k';
    set(ax2,'position',ax1.Position);
    set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
end
%%
print(gcf, 'figures/exFig3.eps', '-depsc', '-vector');