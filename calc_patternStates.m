clear;
addpath(genpath('functions'));

OutputPath = 'data/results/';
DataPath = 'data/fMRI/rfMRI_REST1_LR/';

parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'

[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

File = dir(fullfile(DataPath));
SubjectList = cell(size(File));
for m = 1:length(File)
    SubjectList{m}=File(m).name;
end
SubjectList(1:2)=[];

%%
Start = 101;
T = 300;
DivMap_group = int8(zeros(T*length(SubjectList),nnz(find(V_mask==1))));
CurlMap_group = DivMap_group;
%%
for sub = length(SubjectList):-1:1
    disp(['sub = ' num2str(sub)]);
    L = load([DataPath,SubjectList{sub}, '/VelField_lh_hilbert.mat']);
    R = load([DataPath,SubjectList{sub}, '/VelField_rh_hilbert.mat']);
    
    
    Ux = cat(2,L.VelField.Ux(:,:,Start:Start+T-1),R.VelField.Ux(:,:,Start:Start+T-1));
    Uy = cat(2,L.VelField.Uy(:,:,Start:Start+T-1),R.VelField.Uy(:,:,Start:Start+T-1));
    clearvars L R

    filtSigma = 3;
    filtWidth = round(3*filtSigma);
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    DivMap=[];CurlMap=[];
    for t = T:-1:1

        div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
        div_oneframe(V_mask==0)=nan;
        DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
        
        U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);
        [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
        curl_oneframe(V_mask==0)=nan;
        CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
    end
    clearvars Ux Uy
    
    DivMap(repmat(V_mask,1,1,T)==0)=nan;
    k=2;
    DivMean = mean(DivMap,'all','omitnan');
    DivStd = std(DivMap,[],'all','omitnan');
    DivMap = (DivMap - DivMean)./DivStd;
    DivMap(DivMap<k & DivMap>-k)=0;
    DivMap(DivMap>=k)=1;
    DivMap(DivMap<=-k)=-1;
    DivMap = reshape(DivMap,200*500,[]);
    DivMap = DivMap(V_mask(:)==1,:);
    DivMap_group((sub-1)*T+1:sub*T,:) = int8(DivMap');
    
    CurlMap(repmat(V_mask,1,1,T)==0)=nan;
    CurlMean = mean(CurlMap,'all','omitnan');
    CurlStd = std(CurlMap,[],'all','omitnan');
    CurlMap = (CurlMap - CurlMean)./CurlStd;
    CurlMap(CurlMap<k & CurlMap>-k)=0;
    CurlMap(CurlMap>=k)=1;
    CurlMap(CurlMap<=-k)=-1;
    CurlMap = reshape(CurlMap,200*500,[]);
    CurlMap = CurlMap(V_mask(:)==1,:);
    CurlMap_group(:,(sub-1)*T+1:sub*T) = CurlMap;
end

% save([OutputPath hemisphere '_DivMap_group.mat'],'DivMap_group','-v7.3');
% save([OutputPath hemisphere '_CurlMap_group.mat'],'CurlMap_group','-v7.3');
%%
% Feature = [DivMap_group;CurlMap_group]';
% 
% ref_frame = randperm(size(Feature,1),1);
% for i = size(Feature,1):-1:1
%     if i ~=ref_frame
%         Dist(i) = pdist(Feature([i ref_frame],:));
%     end
% end
% Dist(ref_frame)=[];
% Dist = sort(Dist);
% diffDist = diff(Dist);
% 
% epsilon=15;
% minpts = 1000;
% [idx,corepts] = dbscan(Feature,epsilon,minpts);
% 
% save([OutputPath hemisphere '_DBSCAN_pattern.mat'],'idx','corepts','-v7.3');
%%
 Feature = [DivMap_group;CurlMap_group]';clearvars DivMap_group CurlMap_group

 K = [2 4 8 12];
 for k = 1:length(K)
     [idx,C] = kmeans(Feature,K(k),'Distance','correlation');
     save([OutputPath 'clustering_results/LR_patterns_Kmeans_correlation_K=' num2str(K(k)) '.mat'],'idx','C','-v7.3');
 end