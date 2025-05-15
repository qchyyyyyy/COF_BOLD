clear
clc

addpath(genpath('functions'));
hemisphere = 'LR';

Condition = 'WM';
PE = 'LR';
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks'
[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere(1));
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere(2));

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
%% Load data and task information
load(['data/results/info_task/' Condition '_' PE '_info.mat']);
Block_name = {'0_Back_Body','0_Back_Face',...
    '0_Back_Place','0_Back_Tools',...
    '2_Back_Body','2_Back_Face',...
    '2_Back_Place','2_Back_Tools'};

T = 404;
Events = zeros(T,1);

for bb = 1:8
    idx = ceil(meanStimType(bb,1)-0.5):fix(meanStimType(bb,2)-0.5);
    idx(idx>T)=[];idx(idx<1)=[];
    Events(idx) = meanStimType(bb,3);
end

Events(Events>0 & Events<5)=1;
Events(Events>4)=2;
Events = Events+1;

L = load(['data/results/numPattern_results/trial/lh_numPattern_individual_' Condition '_' PE '_trial_n=105.mat']);
R = load(['data/results/numPattern_results/trial/rh_numPattern_individual_' Condition '_' PE '_trial_n=105.mat']);

SourceMap = cat(2,L.SourceMap,R.SourceMap);
SinkMap = cat(2,L.SinkMap,R.SinkMap);
CCVortexMap = cat(2,L.CCVortexMap,R.CCVortexMap);
CVortexMap = cat(2,L.CVortexMap,R.CVortexMap);

%% temporal correlation
Time = 10:T-9;

Feature = [map2vec(SourceMap(:,:,Time),V_mask);...
    map2vec(SinkMap(:,:,Time),V_mask);...
    map2vec(CCVortexMap(:,:,Time),V_mask);...
    map2vec(CVortexMap(:,:,Time),V_mask)]';

% Feature = [map2vec(SourceMap(:,:,Time),V_mask);...
%     map2vec(SinkMap(:,:,Time),V_mask)]';

% Feature = [map2vec(CCVortexMap(:,:,Time),V_mask);...
%     map2vec(CVortexMap(:,:,Time),V_mask)]';

Feature = zscore(Feature);
Cov = corr(Feature','rows','pairwise');

%% sorted temporal correlation matrix
% Delay = 18;
% Time_Delay = Time - Delay;
% Time_Delay(Time_Delay<1)=[];Time_Delay(Time_Delay>T)=[];
% [sortE,ind] = sort(Events(Time_Delay));
% 
% figure;imagesc(Cov(ind,ind));hold on
% clim([-0.4 0.4]);
% colormap(bluewhitered)
% for i = 1:length(Time)-1
%     if sortE(i)~=sortE(i+1)
%         plot([1 length(Time)],[i+0.5 i+0.5],'LineWidth',3,'Color','k');hold on
%         plot([i+0.5 i+0.5],[1 length(Time)],'LineWidth',3,'Color','k');hold on
%     end
% end
% axis tight;axis equal;
% xlabel('Time (TR)');
% ylabel('Time (TR)');
% % title('Patterns')
% t1=text(18,-12,'fixation');
% t2=text(130,-12,'0 back');
% t3=text(270,-12,'2 back');
%%
rng(0);

[idx,Q] = community_louvain(Cov,1,[],'negative_asym');
K = max(idx);

% K=5;
% Z = linkage(Feature,'single','correlation');
% figure; dendrogram(Z);    
% idx = cluster(Z,K);
%
% Feature = zscore(Feature);
% [idx,C] = kmeans(Feature,K,"Distance","correlation");
% [idx,C] = kmeans(Feature,K,"Distance","correlation",'MaxIter',10000);
% [idx,C] = kmeans(Feature,K,"Distance","correlation",'MaxIter',1000,'Replicates',10,'Options',statset('UseParallel',1));

% decoding accuracy

for Delay = 0:18
    
    Time_Delay = Time-Delay;
    idx0 = idx(length(find(Time_Delay<1))+1:end-length(find(Time_Delay>T)));% cutoff index of clustering results

    Time_Delay(Time_Delay<1)=[];Time_Delay(Time_Delay>T)=[];
    Events_Delay = Events(Time_Delay);% corresponding condition
    
    predict_Events = zeros(size(Events_Delay));

    for k = 1:K
        predict_Events(idx0==k) = mode(Events_Delay(idx0==k));
    end
    decodingAcc(Delay+1)=sum(predict_Events==Events_Delay)./length(Events_Delay);
end
decodingAcc