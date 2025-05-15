clear;
addpath(genpath('functions'));

% delete(gcp('nocreate'));
% partool = parpool('local',35);
% partool.IdleTimeout = Inf;

OutputPath = 'data/results/';
PE = {'LR','RL'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = {'lh','rh'};


for pe = 1
    if pe==1
        Condition = {'MOTOR','LANGUAGE','RELATIONAL','WM','REST1'};% 'WM','RELATIONAL',
    else
        Condition = {'RELATIONAL','WM'};
    end

    for cdn = 1:length(Condition)
        if strcmp(Condition{cdn}(1:end-1),'REST')
            DataPath = ['data/fMRI/rfMRI_REST1_' PE{pe} '/'];
        else
            DataPath = ['data/fMRI/tfMRI_' Condition{cdn} '_' PE{pe} '/'];
        end
        
        
        File = dir(fullfile(DataPath));
        SubjectList = cell(size(File));
        for m = 1:length(File)
            SubjectList{m}=File(m).name;
        end
        SubjectList(1:2)=[];
    for hs = 1:2
        clearvars Source Sink CCVortex CVortex
        switch hemisphere{hs}
            case 'lh'
                [mask, label, oi_x, oi_y, V_mask, mask1] = masklabel(TemplatePath,parcellation,'L');
            case 'rh'
                [mask, label, oi_x, oi_y, V_mask, mask1] = masklabel(TemplatePath,parcellation,'R');
        end

        parfor sub = 1:length(SubjectList)
            disp(['sub = ' num2str(sub)]);
                [Source{sub},Sink{sub},CCVortex{sub},CVortex{sub}] = batch_track_patterns([DataPath,SubjectList{sub}, '/VelField_' hemisphere{hs} '_hilbert.mat'],V_mask);
        end
        
        Source(cellfun(@isempty,Source))=[];
        Sink(cellfun(@isempty,Sink))=[];
        CCVortex(cellfun(@isempty,CCVortex))=[];
        CVortex(cellfun(@isempty,CVortex))=[];
        
        save([OutputPath 'dwelltime_results/' hemisphere{hs} '_PatternDwellTime_group_k=2_' Condition{cdn} '_' PE{pe} '.mat'],...
            'Source','Sink',...
            'CCVortex','CVortex','SubjectList','-v7.3');
    end
    end
end
%% sub function
function [Source,Sink,CCVortex,CVortex] = batch_track_patterns(filename,V_mask)
if exist(filename,'file')==0
    Source=[];Sink=[];CCVortex=[];CVortex=[];
    return
end
load(filename,'VelField');
global T

Start = 16;
T = 200;

Ux = VelField.Ux(:,:,Start:Start+T-1);
Uy = VelField.Uy(:,:,Start:Start+T-1);
clearvars VelField

filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);
    
    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
    
    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
    curl_oneframe(V_mask==0)=nan;
    CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
end

CurlMap(repmat(V_mask,1,1,T)==0)=nan;
DivMap(repmat(V_mask,1,1,T)==0)=nan;

k=2;
DivMean = mean(DivMap,'all','omitnan');
DivStd = std(DivMap,[],'all','omitnan');
DivMap(DivMap<DivMean+k*DivStd & DivMap>DivMean-k*DivStd)=0;
CurlMean = mean(CurlMap,'all','omitnan');
CurlStd = std(CurlMap,[],'all','omitnan');
CurlMap(CurlMap<CurlMean+k*CurlStd & CurlMap>CurlMean-k*CurlStd)=0;
DivMap(DivMap>0)=1;
DivMap(DivMap<0)=-1;
CurlMap(CurlMap>0)=1;
CurlMap(CurlMap<0)=-1;

Source = track_patterns(DivMap,1);% source
Sink = track_patterns(DivMap,-1);% sink
CCVortex = track_patterns(CurlMap,1);%counterclockwise rotation
CVortex = track_patterns(CurlMap,-1);%clockwise rotation
end

function [Source,Sink,CCVortex,CVortex] = batch_track_patterns_surrogate(Path,hemisphere,mask,V_mask)
load([Path, '/TC_surrogate.mat']);
eval(['BrainImg = data.BrainImg_' hemisphere ';']);clearvars data
T0 = size(BrainImg,3)-1;
para_model.tau = 5e-4;
para_model.rho = 1e-8;
BrainImg(isnan(BrainImg))=0;
BrainImg = reshape(BrainImg,[],size(BrainImg,3));
BrainImg = BrainImg(mask(:)==1,:);
BrainImg = transpose(hilbert(BrainImg'));

tmp_BrainImg = zeros([size(mask,1)*size(mask,2) size(BrainImg,2)]);
tmp_BrainImg(mask(:)==1,:) = BrainImg;
clearvars BrainImg
tmp_BrainImg = reshape(tmp_BrainImg,size(mask,1),size(mask,2),[]);
complex_Img{1}=real(tmp_BrainImg);
complex_Img{2} =imag(tmp_BrainImg);
clearvars tmp_BrainImg

global T
% Start = randperm(701,1)+100;
% T = 300;
Start = 101;
T = 1000;

parfor t = 1:T
    [Ux(:,:,t),Uy(:,:,t)] = velocity_field_multidim_ell2({complex_Img{1}(:,:,Start+t-1:Start+t),complex_Img{2}(:,:,Start+t-1:Start+t)},mask,para_model);
end


filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);
    %         div_oneframe= divergence(Ux(:,:,t)./U,Uy(:,:,t)./U);
    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
    
    %     [curl_oneframe,~]= curl(Ux(:,:,t),Uy(:,:,t));
    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
    
    curl_oneframe(V_mask==0)=nan;
    CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
end

CurlMap(repmat(V_mask,1,1,T)==0)=nan;
DivMap(repmat(V_mask,1,1,T)==0)=nan;

k=2;
DivMean = mean(DivMap,'all','omitnan');
DivStd = std(DivMap,[],'all','omitnan');
DivMap(DivMap<DivMean+k*DivStd & DivMap>DivMean-k*DivStd)=0;
CurlMean = mean(CurlMap,'all','omitnan');
CurlStd = std(CurlMap,[],'all','omitnan');
CurlMap(CurlMap<CurlMean+k*CurlStd & CurlMap>CurlMean-k*CurlStd)=0;
DivMap(DivMap>0)=1;
DivMap(DivMap<0)=-1;
CurlMap(CurlMap>0)=1;
CurlMap(CurlMap<0)=-1;

Source = track_patterns(DivMap,1);% source
Sink = track_patterns(DivMap,-1);% sink
CCVortex = track_patterns(CurlMap,1);%counterclockwise rotation
CVortex = track_patterns(CurlMap,-1);%clockwise rotation
end


function DwellTime = track_patterns(Map,value)
global T Cluster
Cluster = zeros(size(Map));
for t = 1:T
    ind = find(Map(:,:,t)==value);
    if length(ind)>=10
        [I,J] = ind2sub(size(Map(:,:,t)),ind);
        Z = linkage([I J],'single');
        K = length(find(Z(:,3)>sqrt(2)))+1;
        Cluster_ind = cluster(Z,'MaxClust',K);
        Cluster_map = zeros(size(Cluster(:,:,t)));
        Cluster_map(ind) = Cluster_ind;
        Cluster(:,:,t) = Cluster_map;
    end
end
clearvars Map Cluster_map Cluster_ind Z K I J

global visited
visited = cell([T 1]);
for t = 1:T
    K = max(Cluster(:,:,t),[],'all');
    visited{t} = zeros([K 1]);
end
Patterns.ind=[];
DwellTime = [];
cnt = 1;
for t = 1:T-1
    K = max(Cluster(:,:,t),[],'all');
    for cc = 1:K
        if ~visited{t}(cc)
            ind_t = find(Cluster(:,:,t)==cc);
            Ind{1} = ind_t;
            Patterns.ind{cnt} = match_cluster(Ind,cc,t);
            DwellTime(cnt) = length(Patterns.ind{cnt});
%             Patterns.time{cnt} = [t t+length(Patterns.ind{cnt})-1];
            cnt = cnt + 1;
        end
    end
end
clearvars Cluster
end % end function

function Ind = match_cluster(Ind,cluster,t)
global T Cluster visited
visited{t}(cluster)=1;
if t == T
    return;
end
for nextcluster = 1:max(Cluster(:,:,t+1),[],'all')
    if ~visited{t+1}(nextcluster)
        ind_nextcluster = find(Cluster(:,:,t+1)==nextcluster);
        if ~isempty(intersect(Ind{end},ind_nextcluster))
            Ind{length(Ind)+1} = ind_nextcluster;
            Ind = match_cluster(Ind,nextcluster,t+1);
            break;
        end
    end
end

end
