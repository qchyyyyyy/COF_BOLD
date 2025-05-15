clear;clc

addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

hcpTemplatePath = 'data/templates/HCP_S1200_Atlas_Z4_pkXDZ/';

flatSurf_L = gifti([hcpTemplatePath 'S1200.L.flat.32k_fs_LR.surf.gii']);
flatSurf_L.vertices(:,1) = flatSurf_L.vertices(:,1) - 250;
flatSurf_R = gifti([hcpTemplatePath 'S1200.R.flat.32k_fs_LR.surf.gii']);
flatSurf_R.vertices(:,1) = flatSurf_R.vertices(:,1) + 250;
flatVertices = [flatSurf_L.vertices;flatSurf_R.vertices];

[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];mask(:,251)=[];
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

load('data/BAT/cifti_database_data.mat');
for termId = length(term_surface_map_res_all):-1:1
    tmp_name = split(dir_path(termId).name,'_pAgF');
    term_name{termId,1} = tmp_name{1};
end

Condition = {'WM','RELATIONAL','LANGUAGE','MOTOR'};
term_taskspecfic = {'encoding_retrieval','matching_task','listening','motor_imagery'};
PE = 'LR';

pattern_names = {'Diverging', 'Converging', 'CC-Spiral', 'C-Spiral'};
pattern_names_old = {'Source', 'Sink', 'CCVortex', 'CVortex'};

% figure('color','w')
set(gcf,'Units','centimeters');
set(gcf,'Position',[1,1,18,18.5]);

t0 = tiledlayout(2,2,"TileSpacing","tight","Padding",'tight');
set(gca,'FontName','Arial','FontUnits','points','FontSize',6,'Visible','off');

for cdn = 1:4
    Task = []; REST = [];
    load(['data/results/' Condition{cdn} '_' PE '_SubList.mat']);
    Task.SubList = SubjectList;
    load('data/results/REST1_LR_SubList.mat');
    REST.SubList = SubjectList;
    SubjectList = intersect(Task.SubList,REST.SubList);
    % voxel level

    Subind = find(ismember(Task.SubList,SubjectList));
    L = load(['data/results/numPattern_results/voxel/' hemisphere{1} '_numPattern_individual_' Condition{cdn} '_' PE '_voxel.mat']);
    R = load(['data/results/numPattern_results/voxel/' hemisphere{2} '_numPattern_individual_' Condition{cdn} '_' PE '_voxel.mat']);
    Task.SourceMap = cat(2,L.SourceMap(:,:,Subind),R.SourceMap(:,:,Subind));
    Task.SinkMap = cat(2,L.SinkMap(:,:,Subind),R.SinkMap(:,:,Subind));
    Task.CCVortexMap = cat(2,L.CCVortexMap(:,:,Subind),R.CCVortexMap(:,:,Subind));
    Task.CVortexMap = cat(2,L.CVortexMap(:,:,Subind),R.CVortexMap(:,:,Subind));


    Subind = find(ismember(REST.SubList,SubjectList));
    L= load(['data/results/numPattern_results/voxel/' hemisphere{1} '_numPattern_individual_REST1_LR_voxel.mat']);
    R = load(['data/results/numPattern_results/voxel/' hemisphere{2} '_numPattern_individual_REST1_LR_voxel.mat']);

    REST.SourceMap = cat(2,L.SourceMap(:,:,Subind),R.SourceMap(:,:,Subind));
    REST.SinkMap = cat(2,L.SinkMap(:,:,Subind),R.SinkMap(:,:,Subind));
    REST.CCVortexMap = cat(2,L.CCVortexMap(:,:,Subind),R.CCVortexMap(:,:,Subind));
    REST.CVortexMap = cat(2,L.CVortexMap(:,:,Subind),R.CVortexMap(:,:,Subind));

    % pair t-test

    nexttile(t0)
    set(gca,'Visible','off');

    t1 = tiledlayout(t0,3,1,"TileSpacing","none","Padding",'none');
    t1.Layout.Tile = cdn;
    
    nexttile(t1)
    set(gca,'Visible','off');

    t2 = tiledlayout(t1,2,2,"TileSpacing","none","Padding","none");
    t2.Layout.Tile = 1;
    t2.Layout.TileSpan = [2, 1];

    Tmap = [];
    Order = [1 3 2 4];
    for i_pattern = 1:length(pattern_names)
        eval(['X = map2vec(Task.' pattern_names_old{i_pattern} 'Map,V_mask);']);
        eval(['Y = map2vec(REST.' pattern_names_old{i_pattern} 'Map,V_mask);']);

        [~,p,~,stats] = ttest(X,Y,"Dim",2);
        t = stats.tstat;
        
        p = FDR_correction(p);
        p = vec2map(p,V_mask);t = vec2map(t,V_mask);
        p(V_mask==0)=nan;t(V_mask==0)=nan;
        Tmap(:,:,i_pattern) = t;

        nexttile(t2);

        imagesc_brainimg(t,p<0.05,0);clim([- 10 10]);
        title(pattern_names{i_pattern});
        hold on
        scatter_boundary(oi_x,oi_y,1);
        
        xlim([10,490]);
        set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
    end

    colormap(bluewhitered);

    % functional enrichment analysis
    for i_pattern = 1:length(pattern_names)
        eval(['contrast_' pattern_names_old{i_pattern} ' = mean(Task.' pattern_names_old{i_pattern} 'Map,3)-mean(REST.' pattern_names_old{i_pattern} 'Map,3);']);
        eval(['contrast_patterns(:,:,' num2str(i_pattern) ') = mean(Task.' pattern_names_old{i_pattern} 'Map,3)-mean(REST.' pattern_names_old{i_pattern} 'Map,3);']);
    end

    % term_name = [];T = [];
    % for termId = length(term_surface_map_res_all):-1:1
    %     cData = term_surface_map_res_all{termId}(:,2);
    %     idx = term_surface_map_res_all{termId}(:,1);
    %     TC = zeros(size(flatVertices(:,1)));
    %     TC(idx) = cData;
    % 
    %     Map = tc2mesh(TC,flatVertices(:,1),flatVertices(:,2),mask);
    %     % Map(Map>0)=1;Map(Map<0)=1;
    %     V_Map = conv2(Map,ones(2,2)/4,'valid');
    % 
    %     % figure;
    %     idx_rela = find(V_Map>0 & V_mask==1);
    %     idx_irrela = setdiff(find(V_mask==1),idx_rela);
    %     for i_pattern = length(pattern_names):-1:1
    %         t = Tmap(:,:,i_pattern);
    %         [~,P(termId,i_pattern),~,stats] = ttest2(t(idx_rela),t(idx_irrela),'Vartype','unequal');
    %         T(termId,i_pattern) = stats.tstat;
    %     end
    % end

    termId = find(contains(term_name,term_taskspecfic{cdn}));
    cData = term_surface_map_res_all{termId}(:,2);
    idx = term_surface_map_res_all{termId}(:,1);
    TC = zeros(size(flatVertices(:,1)));
    TC(idx) = cData;

    Map = tc2mesh(TC,flatVertices(:,1),flatVertices(:,2),mask);

    V_Map = conv2(Map,ones(2,2)/4,'valid');


    idx_rela = find(V_Map>0 & V_mask==1);
    idx_irrela = setdiff(find(V_mask==1),idx_rela);

    nexttile(t1)
    set(gca,'Visible','off');

    t2 = tiledlayout(t1,1,4,"TileSpacing","compact","Padding","compact");
    t2.Layout.Tile = 3;
    for i_pattern = 1:length(pattern_names)
        t = Tmap(:,:,i_pattern);

        nexttile(t2)
        g1 = repmat({'Task-relevant'},length(idx_rela),1);
        g2 = repmat({'Task-irrelevant'},length(idx_irrela),1);
        g = [g1; g2];

        boxplot([t(idx_rela);t(idx_irrela)],[g1; g2],'Symbol', '');
        [~,PPP{cdn,i_pattern},~,SSS{cdn,i_pattern}] = ttest2(t(idx_rela),t(idx_irrela),'Vartype','unequal');
        ylim([-15 20]);
        title(pattern_names{i_pattern});
        ylabel('t-value');

        set(gca,'FontName','Arial','FontUnits','points','FontSize',6);
    end

end
%%
% print(gcf, 'figures/exFig8.eps', '-depsc', '-vector');