clear;clc
addpath(genpath('functions'))

dataPath = 'data/fMRI/rfMRI_REST1_LR/';
% dataPath = 'data/fMRI/tfMRI_WM_LR/';
% dataPath = 'data/fMRI/tfMRI_WM_RL/';
% dataPath = 'data/fMRI/tfMRI_RELATIONAL_LR/';
% dataPath = 'data/fMRI/tfMRI_RELATIONAL_RL/';
% dataPath = 'data/fMRI/tfMRI_LANGUAGE_LR/';
% dataPath = 'data/fMRI/tfMRI_MOTOR_LR/';

File = dir(fullfile(dataPath));
SubjectList = cell(size(File));
for m = 1:length(File)
   SubjectList{m}=File(m).name;
end
SubjectList(1:3)=[];
%% generation and preprocessing of the surrogate data
surrogateMethod = '2D';
for sub = 1:length(SubjectList)
    disp(['sub = ' num2str(sub)]);
    preproc_2D_surrogate([dataPath,SubjectList{sub},'/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'],...
        [dataPath,SubjectList{sub},'/'],surrogateMethod);
end

%% VVFs estimation via COF
load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_mask_right');
tau = 5e-4;rho = 1e-8;
% tau = 1e-9;rho = 1;
sigma = 8;

parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(parcellation,'R');

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
% clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R


Start = 201;
T = 100;
sub = 1;
load([dataPath,SubjectList{sub},'/TC_surrogate_2D.mat']);
surrogateData = data.BrainImg_lh(:,:,Start:Start+T);

[Ux_s,Uy_s] = extract_velocity_field_hilbert(surrogateData,cii_mask_left,tau,rho,sigma);

load([dataPath,SubjectList{sub},'/TC_mesh_r_f_zscore.mat']);
realData = data.BrainImg_lh(:,:,Start:Start+T);

[Ux_r,Uy_r] = extract_velocity_field_hilbert(realData,cii_mask_left,tau,rho,sigma);


% sigma = 1;
% [Ux_s_unfilt,Uy_s_unfilt] = extract_velocity_field_hilbert(surrogateData,cii_mask_left,tau,rho,sigma);
% [Ux_r_unfilt,Uy_r_unfilt] = extract_velocity_field_hilbert(realData,cii_mask_left,tau,rho,sigma);
% 
% [r,p] = corr([Ux_s_unfilt(repmat(V_mask_L,1,1,T)==1);Uy_s_unfilt(repmat(V_mask_L,1,1,T)==1)],...
%     [Ux_s(repmat(V_mask_L,1,1,T)==1);Uy_s(repmat(V_mask_L,1,1,T)==1)])
% 
% [r,p] = corr([Ux_r_unfilt(repmat(V_mask_L,1,1,T)==1);Uy_r_unfilt(repmat(V_mask_L,1,1,T)==1)],...
%     [Ux_r(repmat(V_mask_L,1,1,T)==1);Uy_r(repmat(V_mask_L,1,1,T)==1)])
%%
filtSigma = sigma;
filtWidth = ceil(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

F(T) = struct('cdata',[],'colormap',[]);
figure('color','w');
for t = 1:T
    subplot(1,2,1);
    imagesc_brainimg(realData(:,:,t),mask_L,0);colorbar();clim([-3 3]);hold on
    % imagesc_brainimg(nanconv(realData(:,:,t),imageFilter),mask_L,0);colorbar();clim([-3 3]);hold on
    % quiver_dsr(Ux_r_unfilt(:,:,t),Uy_r_unfilt(:,:,t),2,0);hold off
    qq = quiver_dsr(Ux_r(:,:,t),Uy_r(:,:,t),2,0);qq.Color='k';hold off
    subplot(1,2,2);
    imagesc_brainimg(surrogateData(:,:,t),mask_L,0);colorbar();clim([-3 3]);hold on
    % imagesc_brainimg(nanconv(surrogateData(:,:,t),imageFilter),mask_L,0);colorbar();clim([-3 3]);hold on
    % quiver_dsr(Ux_s_unfilt(:,:,t),Uy_s_unfilt(:,:,t),2,0);hold off
    qq = quiver_dsr(Ux_s(:,:,t),Uy_s(:,:,t),2,0);qq.Color='k';hold off
    pause(0.1);
end
%%
do_zscore = 0;
filtSigma=3;
[DivMap_r, CurlMap_r] = calc_DivCurl(Ux_r(:,:,11:90),Uy_r(:,:,11:90),V_mask_L,filtSigma,do_zscore);
[DivMap_s, CurlMap_s] = calc_DivCurl(Ux_s(:,:,11:90),Uy_s(:,:,11:90),V_mask_L,filtSigma,do_zscore);

% MeanDivMap_s = mean(DivMap_s,'all','omitnan')
% StdDivMap_s = std(DivMap_s,[],'all','omitnan')

% %%
% [r,p] = corr([Ux_r_unfilt(5e5+DivMap_r>2);Uy_r_unfilt(5e5+DivMap_r>2)],...
%     [Ux_r(5e5+DivMap_r>2);Uy_r(5e5+DivMap_r>2)])
% 
% [r,p] = corr([Ux_s_unfilt(5e5+DivMap_s>2);Uy_s_unfilt(5e5+DivMap_s>2)],...
%     [Ux_s(5e5+DivMap_s>2);Uy_s(5e5+DivMap_s>2)])
% 
% [r,p] = corr([Ux_r_unfilt(5e5+CurlMap_r>2);Uy_r_unfilt(5e5+CurlMap_r>2)],...
%     [Ux_r(5e5+CurlMap_r>2);Uy_r(5e5+CurlMap_r>2)])
% 
% [r,p] = corr([Ux_s_unfilt(5e5+CurlMap_s>2);Uy_s_unfilt(5e5+CurlMap_s>2)],...
%     [Ux_s(5e5+CurlMap_s>2);Uy_s(5e5+CurlMap_s>2)])
%% figure 4a
[h, p, ci, stats] = vartest2(DivMap_r(repmat(V_mask_L,1,1,80)==1), DivMap_s(repmat(V_mask_L,1,1,80)==1))
[h, p, ci, stats] = vartest2(CurlMap_r(repmat(V_mask_L,1,1,80)==1), CurlMap_s(repmat(V_mask_L,1,1,80)==1))


alphaValue = 0.3;
% figure;
clf
subplot(2,1,1);
hh = histogram(DivMap_r(repmat(V_mask_L,1,1,80)==1),-2:0.05:2,'Normalization','probability');hold on
alpha(hh,alphaValue);
hh = histogram(DivMap_s(repmat(V_mask_L,1,1,80)==1),-2:0.05:2,'Normalization','probability');
alpha(hh,alphaValue);

legend('original','surrogate','FontSize',12);
xlim([-0.8 0.8])
ax=gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ylabel('pdf','FontSize',14);
xlabel('divergence','FontSize',14);

subplot(2,1,2);
hh = histogram(CurlMap_r(repmat(V_mask_L,1,1,80)==1),-2:0.01:2,'Normalization','probability');hold on
alpha(hh,alphaValue);
hh = histogram(CurlMap_s(repmat(V_mask_L,1,1,80)==1),-2:0.01:2,'Normalization','probability');
alpha(hh,alphaValue);
legend('original','surrogate','FontSize',12);
xlim([-0.5 0.5])
ax=gca;
ax.LineWidth = 1.2;
ax.FontSize = 12;
ylabel('pdf','FontSize',14);
xlabel('normalized curl','FontSize',14);

set(gcf,'Position',[543.5,114,456.5,584]);
%%
F(T) = struct('cdata',[],'colormap',[]);
figure('color','w');
for t = 1:T
    subplot(1,2,1);
    imagesc_brainimg(DivMap_r(:,:,t),V_mask_L,0);colorbar();clim([-3 3]);hold on
    quiver_dsr(Ux_r(:,:,t+10),Uy_r(:,:,t+10),2,0);hold off
    subplot(1,2,2);
    imagesc_brainimg(DivMap_s(:,:,t),V_mask_L,0);colorbar();clim([-3 3]);hold on
    quiver_dsr(Ux_s(:,:,t+10),Uy_s(:,:,t+10),2,0);hold off
    pause(1);
end