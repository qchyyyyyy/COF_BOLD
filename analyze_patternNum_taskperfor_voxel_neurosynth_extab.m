clear;clc
addpath(genpath('functions'));
hemisphere = {'lh','rh'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

Condition = 'WM';
PE = 'RL';
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];mask(:,251)=[];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
load(['data/results/info_task/' Condition '_' PE '_info.mat']);

%%
hcpTemplatePath = 'data/templates/HCP_S1200_Atlas_Z4_pkXDZ/';

flatSurf_L = gifti([hcpTemplatePath 'S1200.L.flat.32k_fs_LR.surf.gii']);
flatSurf_L.vertices(:,1) = flatSurf_L.vertices(:,1) - 250;
flatSurf_R = gifti([hcpTemplatePath 'S1200.R.flat.32k_fs_LR.surf.gii']);
flatSurf_R.vertices(:,1) = flatSurf_R.vertices(:,1) + 250;
flatVertices = [flatSurf_L.vertices;flatSurf_R.vertices];
clearvars flatSurf_L flatSurf_R hcpTemplatePath
%%
Num_L = load(['data/results/numPattern_results_new/index/' hemisphere{1} '_numPattern_individual_' Condition '_' PE '_index.mat']);
Num_R = load(['data/results/numPattern_results_new/index/' hemisphere{2} '_numPattern_individual_' Condition '_' PE '_index.mat']);

Num_R.Source = cellfun(@(x) x+200*250,Num_R.Source,'UniformOutput',false);
Num_R.Sink = cellfun(@(x) x+200*250,Num_R.Sink,'UniformOutput',false);
Num_R.CCVortex = cellfun(@(x) x+200*250,Num_R.CCVortex,'UniformOutput',false);
Num_R.CVortex = cellfun(@(x) x+200*250,Num_R.CVortex,'UniformOutput',false);

load('data/BAT/cifti_database_data.mat','term_surface_map_res_all','dir_path');

%%
for termId = length(term_surface_map_res_all):-1:1
    disp(['termId = ' num2str(termId)]);
    cData = term_surface_map_res_all{termId}(:,2);
    idx = term_surface_map_res_all{termId}(:,1);
    TC = zeros(size(flatVertices(:,1)));
    TC(idx) = cData;

    Map = tc2mesh(TC,flatVertices(:,1),flatVertices(:,2),mask);

    V_Map = conv2(Map,ones(2,2)/4,'valid');


    idx_rela = find(V_Map>0 & V_mask==1);

    % idx_irrela = setdiff(find(V_mask==1),idx_rela);
    % clearvars termId cData idx TC Map V_Map flatVertices
    

    NumSource = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.Source, Num_R.Source, 'UniformOutput',false);
    NumSource = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumSource, 'UniformOutput', false));
    NumSink = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.Sink, Num_R.Sink, 'UniformOutput',false);
    NumSink = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumSink, 'UniformOutput', false));
    NumCCVortex = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.CCVortex, Num_R.CCVortex, 'UniformOutput',false);
    NumCCVortex = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumCCVortex, 'UniformOutput', false));
    NumCVortex = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.CVortex, Num_R.CVortex, 'UniformOutput',false);
    NumCVortex = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumCVortex, 'UniformOutput', false));
    % clearvars Num_L Num_R
    
    Subind=1:size(NumSource,1);
    T = 5:size(NumSource,2)-4;
    T(T>size(NumSource,2))=[];
    NumSource = NumSource(Subind,T);
    NumSink = NumSink(Subind,T);
    NumCCVortex = NumCCVortex(Subind,T);
    NumCVortex = NumCVortex(Subind,T);

    Performance = Acc;
    Performance = Performance(Subind,~isnan(Performance(1,:)));
    Threshold = 0;

    
    % YYLim = [600 1400];
    YYLim = [0.3 1.2];
    patternNames = {'Source','Sink','CCVortex','CVortex'};

    % figure;
    % tt = tiledlayout(2,5,"TileSpacing","compact","Padding",'compact');
    % xlabel(tt,'Number of Patterns','FontSize',15);
    % ylabel(tt,'Accuracy','FontSize',15);
    % title(tt,'Working Memory','FontSize',18);
    % title(tt,'Relational Processing','FontSize',18);

    for i_pattern = 1:4
        nexttile
        eval(['x=mean(Num' patternNames{i_pattern} ',2);'])
        [x,idx] = sort(x);
        y=mean(Performance,2);y = y(idx);
        x = x(y>0)';y = y(y>0)';

        [R(termId,i_pattern),P(termId,i_pattern)] = corr(x',y');

        % [b,stats] = robustfit(x,y);p = flip(b);t(i_pattern,1) = stats.t(2);
        % y1=polyval(p,x);
        % yfit=polyconf(p,x);
        % alpha = 0.05;
        % t_value = tinv(1 - alpha/2, stats.dfe);
        % se_pred = sqrt(stats.s^2 * (1 + sum(([ones(size(x')) x'] / stats.R).^2, 2)));
        % dy = t_value * se_pred';
        %
        % h1=fill([x,fliplr(x)],[yfit-dy,fliplr(yfit+dy)],[0.8706 0.9216 0.9804]);
        % hold on
        % h2=plot(x,y1,'b','linewidth',2);
        % hold on
        % h3=plot(x,y,'r.','markersize',15);
        %
        % xlim([min([x,fliplr(x)]) max([x,fliplr(x)])]);ylim(YYLim);
        %
        % title(patternNames{i_pattern});
        % if i_pattern==1;ylabel('Run 1','FontSize',10);end
        %
        % % ylabel('Surrogate','FontSize',10);
        % text(0.95*min([x,fliplr(x)])+0.05*max([x,fliplr(x)]),YYLim(2)-0.05,{['P = ' num2str(stats.p(2))]});
    end

    nexttile
    x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);[x,idx] = sort(x);
    y=mean(Performance,2);y = y(idx);
    x = x(y>0)';y = y(y>0)';% x = x(y~=0)';y = y(y~=0)';

    [R(termId,5),P(termId,5)] = corr(x',y');

    % [b,stats] = robustfit(x,y);p = flip(b);t(5,1) = stats.t(2);
    % y1=polyval(p,x);
    % yfit=polyconf(p,x);
    % alpha = 0.05;
    % t_value = tinv(1 - alpha/2, stats.dfe);
    % se_pred = sqrt(stats.s^2 * (1 + sum(([ones(size(x')) x'] / stats.R).^2, 2)));
    % dy = t_value * se_pred';
    % h1=fill([x,fliplr(x)],[yfit-dy,fliplr(yfit+dy)],[0.8706 0.9216 0.9804]);
    % hold on
    % h2=plot(x,y1,'b','linewidth',2);
    % hold on
    % h3=plot(x,y,'r.','markersize',15);
    % xlim([min([x,fliplr(x)]) max([x,fliplr(x)])]);ylim(YYLim);
    % title('All');
    % text(0.95*min([x,fliplr(x)])+0.05*max([x,fliplr(x)]),YYLim(2)-0.05,{['P = ' num2str(stats.p(2))]});
    %
    % set(gcf,'Position',[2.5,27,1408,331*2]);

end

[~,idx] = sort(mean(R(:,[1 2 5]),2),'descend');
for termId = 1:length(term_surface_map_res_all)
    tmp_name = split(dir_path(termId).name,'_pAgF');
    term_name{termId,1} = tmp_name{1};
end
term_name = term_name(idx);
R = R(idx,:);
P = P(idx,:);
%%
patternNames = {'Source','Sink','CCVortex','CVortex'};
idx_rela = find(V_mask==1);
NumSource = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.Source, Num_R.Source, 'UniformOutput',false);
NumSource = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumSource, 'UniformOutput', false));
NumSink = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.Sink, Num_R.Sink, 'UniformOutput',false);
NumSink = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumSink, 'UniformOutput', false));
NumCCVortex = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.CCVortex, Num_R.CCVortex, 'UniformOutput',false);
NumCCVortex = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumCCVortex, 'UniformOutput', false));
NumCVortex = cellfun(@(x, y) cat(1, x(:), y(:)), Num_L.CVortex, Num_R.CVortex, 'UniformOutput',false);
NumCVortex = cell2mat(cellfun(@(z) nnz(z(ismember(z, idx_rela))), NumCVortex, 'UniformOutput', false));

Subind=1:size(NumSource,1);
T = 5:size(NumSource,2)-4;
T(T>size(NumSource,2))=[];
NumSource = NumSource(Subind,T);
NumSink = NumSink(Subind,T);
NumCCVortex = NumCCVortex(Subind,T);
NumCVortex = NumCVortex(Subind,T);

Performance = Acc;
Performance = Performance(Subind,~isnan(Performance(1,:)));
Threshold = 0;
%%
for i_pattern = 1:4
    eval(['x=mean(Num' patternNames{i_pattern} ',2);'])
    [x,idx] = sort(x);
    y=mean(Performance,2);y = y(idx);
    x = x(y>0)';y = y(y>0)';
    [RR(i_pattern),PP(i_pattern)] = corr(x',y');
end

nexttile
x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);[x,idx] = sort(x);
y=mean(Performance,2);y = y(idx);
x = x(y>0)';y = y(y>0)';% x = x(y~=0)';y = y(y~=0)';

[RR(5),PP(5)] = corr(x',y');
%%
% PE = 'LR';
%
% if length(hemisphere)==2
%     Num_L = load([hemisphere{1} '_numPattern_individual_' Condition '_' PE '.mat']);
%     Num_R = load([hemisphere{2} '_numPattern_individual_' Condition '_' PE '.mat']);
%     NumSource = Num_L.NumSource + Num_R.NumSource;
%     NumSink = Num_L.NumSink + Num_R.NumSink;
%     NumCCVortex = Num_L.NumCCVortex + Num_R.NumCCVortex;
%     NumCVortex = Num_L.NumCVortex + Num_R.NumCVortex;
% else
%     load([hemisphere{1} '_numPattern_individual_' Condition '_' PE '.mat']);
% end
% Subind=1:size(NumSource,1);
% NumSource = NumSource(Subind,5:end-4);
% NumSink = NumSink(Subind,5:end-4);
% NumCCVortex = NumCCVortex(Subind,5:end-4);
% NumCVortex = NumCVortex(Subind,5:end-4);
%
% load([Condition '_' PE '_info.mat']);
%
% RT = RT(Subind,~isnan(RT(1,:)));
%
% %%
%
% for i_pattern = 1:4
%     nexttile
%     eval(['x=mean(Num' patternNames{i_pattern} ',2);'])
%     [x,idx] = sort(x);
%     y=mean(RT,2);y = y(idx);
%     x = x(y>0)';y = y(y>0)';
%     [b,stats] = robustfit(x,y);p = flip(b);t(i_pattern,2) = stats.t(2);
%     y1=polyval(p,x);
%     yfit=polyconf(p,x);
%     alpha = 0.05;
%     t_value = tinv(1 - alpha/2, stats.dfe);
%     se_pred = sqrt(stats.s^2 * (1 + sum(([ones(size(x')) x'] / stats.R).^2, 2)));
%     dy = t_value * se_pred';
%
%     h1=fill([x,fliplr(x)],[yfit-dy,fliplr(yfit+dy)],[0.8706 0.9216 0.9804]);
%     hold on
%     h2=plot(x,y1,'b','linewidth',2);
%     hold on
%     h3=plot(x,y,'r.','markersize',15);
%
%     xlim([min([x,fliplr(x)]) max([x,fliplr(x)])]);ylim(YYLim);
%
%     title(patternNames{i_pattern});
%     if i_pattern==1;ylabel('Run 2','FontSize',10);end
%
%     % ylabel('Surrogate','FontSize',10);
%     text(0.95*min([x,fliplr(x)])+0.05*max([x,fliplr(x)]),YYLim(2)-0.05,{['P = ' num2str(stats.p(2))]});
% end
%
% nexttile
% x=mean(NumSource,2)+mean(NumSink,2)+mean(NumCCVortex,2)+mean(NumCVortex,2);[x,idx] = sort(x);
% y=mean(RT,2);y = y(idx);
% x = x(y>0)';y = y(y>0)';% x = x(y~=0)';y = y(y~=0)';
% [b,stats] = robustfit(x,y);p = flip(b);t(5,2) = stats.t(2)
% y1=polyval(p,x);
% yfit=polyconf(p,x);
% alpha = 0.05;
% t_value = tinv(1 - alpha/2, stats.dfe);
% se_pred = sqrt(stats.s^2 * (1 + sum(([ones(size(x')) x'] / stats.R).^2, 2)));
% dy = t_value * se_pred';
% h1=fill([x,fliplr(x)],[yfit-dy,fliplr(yfit+dy)],[0.8706 0.9216 0.9804]);
% hold on
% h2=plot(x,y1,'b','linewidth',2);
% hold on
% h3=plot(x,y,'r.','markersize',15);
% xlim([min([x,fliplr(x)]) max([x,fliplr(x)])]);ylim(YYLim);
% title('All');
% text(0.95*min([x,fliplr(x)])+0.05*max([x,fliplr(x)]),YYLim(2)-0.05,{['P = ' num2str(stats.p(2))]});
%
% set(gcf,'Position',[2.5,27,1408,331*2]);