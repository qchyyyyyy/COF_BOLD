clear;

addpath(genpath('functions'));
DataPath = 'data/fMRI/rfMRI_REST1_LR/';
OutputPath = 'data/results/decomposition_results/';

parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(TemplatePath,parcellation,hemisphere(1));
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(TemplatePath,parcellation,hemisphere(2));

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

dsr = 3;
mask_s = zeros(size(mask));
mask_s(1:dsr:end,1:dsr:end) = 1;mask_s(mask==0)=0;
%%
Nnodes = nnz(find(mask_s==1));
lagMatrix = zeros(Nnodes,Nnodes,length(SubjectList));

delete(gcp('nocreate'));
partool = parpool('local',35);
partool.IdleTimeout = Inf;

parfor sub = 1:length(SubjectList)
    disp(['sub = ' num2str(sub)]);
    lagMatrix(:,:,sub) = generate_lagMatrix(...
        [DataPath,SubjectList{sub},'/'],...
        mask_s);
end
%%
GrouplagMatrix = nanmean(lagMatrix,3);clearvars lagMatrix

save([OutputPath hemisphere '_lagProjection_group_REST.mat'],...
    'GrouplagMatrix','-v7.3');
%%

function LagMatrix = generate_lagMatrix(Path,mask)

load([Path,'TC_mesh_r_f_zsore_s.mat'],'data');
Start = 101;
T = min(1000,size(data.BrainImg_lh,3)-Start-98);

BrainImg = cat(2,data.BrainImg_lh(:,1:end-1,Start:Start+T),data.BrainImg_rh(:,:,Start:Start+T));
clearvars data

BrainImg = map2vec(BrainImg,mask)';

%% the following part is a reconfiguration of previous work (https://github.com/RaichleLab/lag-code)
% Specify data parameters

lag_lim = 4;    % lag limit (in seconds)
lags = -7:7;    % range of TR shifts; max(lags) = round(lag_lim/tr + 1)
TR = 0.72; % sampling interval in seconds
Nframes = size(BrainImg,1);
Nnodes = size(BrainImg,2);

% Do the lagged correlation/covariance computation of TD matrices
BrainImg = zscore(BrainImg);% De-mean time series

% Loop over blocks of contiguous frames
L = (Nframes - abs(lags));
L = reshape(L,1,1,[]);
LagMatrix = zeros(Nnodes,Nnodes);
for i = Nnodes:-1:1
    for j = Nnodes:-1:1
        Cov = lagged_cov(BrainImg(:,i),BrainImg(:,j),max(lags));
        % Normalize pairwise cross-covariance functions based on entire run
        Cov = Cov./L;
        
        % Parabolic interpolation to get peak lag/correlation
        [subj_lags,~] = parabolic_interp(Cov,TR);
        subj_lags(abs(subj_lags) > lag_lim) = nan;% Exclude long lags (generally occur when CCF is flat)
        LagMatrix(i,j) = subj_lags;
    end
end

end

function  r = lagged_cov(Avg1,Avg2,L)
% From: https://github.com/ryraut/lag-code/blob/master/lagged_cov.m
%
% Citation: Mitra, Anish, et al. "Lag structure in resting-state fMRI." Journal of Neurophysiology 111.11 (2014): 2374-2391.
% Raut, Ryan V., et al. "On time delay estimation and sampling error in resting-state fMRI." Neuroimage 194 (2019): 211-227.


L1 = size(Avg1,2);
L2 = size(Avg2,2);
r = single(zeros(L1,L2,2*L+1));

k = 1;
for i = -L:L
    tau = abs(i);
    
    if i >=0
        Avg1_lagged = Avg1(1:end-tau,:);
        Avg2_lagged = Avg2(1+tau:end,:);
    else
        Avg1_lagged = Avg1(1+tau:end,:);
        Avg2_lagged = Avg2(1:end-tau,:);
    end
    
    r(:,:,k) = Avg1_lagged'*Avg2_lagged;
    k = k+1;
end

end

function [peak_lag,peak_cov] = parabolic_interp(lcc,tr)
% From: https://github.com/ryraut/lag-code/blob/master/parabolic_interp.m
%
% Citation: Mitra, Anish, et al. "Lag structure in resting-state fMRI." Journal of Neurophysiology 111.11 (2014): 2374-2391.
% Raut, Ryan V., et al. "On time delay estimation and sampling error in resting-state fMRI." Neuroimage 194 (2019): 211-227.

s = size(lcc);
peak_lag = nan([1,s(1)*s(2)]);
peak_cov = peak_lag;

% linearize
lcc = reshape(lcc,[s(1)*s(2),s(3)])';

% find index of extremum (max or min determined by sign at zero-lag)
[~,I]= max(bsxfun(@times,lcc,sign(lcc((s(3)+1)/2,:))),[],1);

% ensure extremum is not at an endpoint (this would preclude parabolic interpolation)
use = I>1 & I<s(3);
if use==0
    return
end
lcc = lcc(:,use);

% place peaks at center
x0 = I(use) - (s(3)+1)/2;

% set up three-point ccf for interpolation (y1,y2,y3)
i = sub2ind([size(lcc),sum(use)],I(use),1:sum(use));
lcc = [lcc(i-1);lcc(i);lcc(i+1)];

% fit parabola: tau = TR * (y1-y3) / (2*(y1-2y2+y3))
b = (lcc(3,:) - lcc(1,:))/2;
a = (lcc(1,:) + lcc(3,:) - 2*lcc(2,:))/2;
peak_lag(use) =  (-b./(2*a));

% construct parabola to get covariance (y = ax^2 + bx + c)
peak_cov(use) = a.*(peak_lag(use).^2) + b.*peak_lag(use) + lcc(2,:);

% put back TR information
peak_lag(use) = (peak_lag(use) + x0)*tr;

peak_lag = reshape(peak_lag,[s(1) s(2)]);
peak_cov = reshape(peak_cov,[s(1) s(2)]);

end