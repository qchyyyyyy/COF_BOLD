clear;clc
addpath(genpath('functions'));
hemisphere = 'lh';
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks'
switch hemisphere
    case 'lh'
    [mask,label,oi_x,oi_y,V_mask] = masklabel(parcellation,'L');
    case 'rh'
    [mask,label,oi_x,oi_y,V_mask] = masklabel(parcellation,'R');
end
V_label = max(label(1:end-1,1:end-1),max(label(2:end,1:end-1),max(label(1:end-1,2:end),label(2:end,2:end))));
V_label(V_mask==0)=0;

[M,N] = size(mask);
%% Load data
Start = 201;
T = 50;

load('data/fMRI/rfMRI_REST1_LR/120010/VelField_lh_hilbert.mat');
Ux = VelField.Ux(:,:,Start:Start+T-1);
Uy = VelField.Uy(:,:,Start:Start+T-1);
%%
filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);
    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
    curl_oneframe(V_mask==0)=nan;
    CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);

    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
end

CurlMap(repmat(V_mask,1,1,T)==0)=nan;DivMap(repmat(V_mask,1,1,T)==0)=nan;

k=2;
DivMean = mean(DivMap,'all','omitnan');DivStd = std(DivMap,[],'all','omitnan');DivMap = (DivMap-DivMean)./DivStd;
CurlMean = mean(CurlMap,'all','omitnan');CurlStd = std(CurlMap,[],'all','omitnan');CurlMap = (CurlMap-CurlMean)./CurlStd;
%%

F(T) = struct('cdata',[],'colormap',[]);
[vectorImg, ~, alphaValue] = imread('data/figures/colormap_patterns.png');

figure('color','k');
for t = 1:T
    clf;
    
%     ax1=axes;
%     q = quiver_dsr(Ux(:,:,t),Uy(:,:,t),2,0);
%     set(q,'Color','w');
%     ax1.Visible='off';
%     axis([0 250 0 200]);axis equal;axis tight;

    ax2 = axes;
    imagesc_brainimg(CurlMap(:,:,t),~isnan(CurlMap(:,:,t)),0,[-3 3]);%,[CurlMean-k*CurlStd CurlMean+k*CurlStd]
    alpha(abs(CurlMap(:,:,t))/3);
    colormap(ax2,winter);
    ax2.Visible='off';
    axis([0 250 0 200]);axis equal;axis tight;

    ax3 = axes;
    imagesc_brainimg(DivMap(:,:,t),~isnan(DivMap(:,:,t)),0,[-3 3]);%,[DivMean-k*DivStd DivMean+k*DivStd]
    alpha(abs(DivMap(:,:,t))/3);
    colormap(ax3,autumn);
    ax3.Visible='off';
    axis([0 250 0 200]);axis equal;axis tight;

    hold on
    q = streamslice(Ux(:,:,t),Uy(:,:,t),16);

    set(q,'Color','w');
    hold off
    
    ax4 = axes;
    set(ax4,'Position',[0.6 0 0.4 0.2]);
    K=10;
    xLim = get(gca, 'XLim');
    imshow(vectorImg);
    axis off
    
    drawnow
    F(t) = getframe(gcf);
end
createmovie(F,'data/videos/SupplemantaryVideo2',1./0.72);