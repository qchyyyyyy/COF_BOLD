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
T = 100;

load('data/fMRI/rfMRI_REST1_LR/120010/TC_mesh_r_f_zscore.mat');
switch hemisphere
    case 'lh'
    BrainImg = data.BrainImg_lh(:,:,Start:Start+T);
    case 'rh'
    BrainImg = data.BrainImg_rh(:,:,Start:Start+T);
end

BrainImg(repmat(mask,1,1,T+1)==0)=nan;
MeanImg = mean(BrainImg,'all','omitnan');
StdImg = std(BrainImg,[],'all','omitnan');

analyticBrainImg = map2vec(BrainImg,mask);%convert size of BrainImg from M*N*T to voxel*T
analyticBrainImg = hilbert(analyticBrainImg').';% hilbert transformation
analyticBrainImg = vec2map(analyticBrainImg,mask);% convert size of BrainImg from voxel*T to M*N*T

load('data/fMRI/rfMRI_REST1_LR/120010/VelField_lh_hilbert.mat');
Ux = VelField.Ux(:,:,Start:Start+T-1);
Uy = VelField.Uy(:,:,Start:Start+T-1);

%% extended data figure 1
figure;
% subplot(2,1,1);% figure 2a-1
pp = patch([1:101 NaN],[squeeze(real(analyticBrainImg(100,100,:)))' NaN],...
    [squeeze(angle(analyticBrainImg(100,100,:)))' NaN],'Marker','none','EdgeColor','interp','MarkerFaceColor','flat');
pp.LineWidth=2.5;
colormap(hsv)
clim([-pi pi])
xlim([1 101]);
ylim([-3 3]);
xlabel('Time (TR)','FontSize',16);
ylabel('fMRI signal','FontSize',16);
ax=gca;ax.LineWidth=1;
set(gcf,'Position',[452.5,206.5,396.5,394.5/2]);
set(gca,'FontSize',15);
exportgraphics(gcf,'figures/exFig1_1.eps','Colorspace','rgb','Resolution',600);

%%
figure;% subplot(2,1,2);% figure 2a-2
pp = patch([1:101 NaN],[squeeze(real(analyticBrainImg(100,100,:)))' NaN],...
    [squeeze(imag(analyticBrainImg(100,100,:)))' NaN],...
    [squeeze(angle(analyticBrainImg(100,100,:)))' NaN],'Marker','none','EdgeColor','interp','MarkerFaceColor','flat');
pp.LineWidth=2.5;
colormap(hsv)
clim([-pi pi])
xlim([1 101]);
%xl = xlabel('Time (TR)','FontSize',16);xl.Rotation=5;
%yl = ylabel('Real','FontSize',16);yl.Rotation=-40;
%zlabel('Imaginary','FontSize',16);
ax=gca;ax.LineWidth=1;
ax.XTickLabel=[];
ax.YTickLabel=[];
ax.ZTickLabel=[];

view(-10,30);
set(gca,'FontName','Arial');
set(gcf,'Position',[452.5,206.5,396.5,394.5/2]);
set(gca,'FontSize',15);

exportgraphics(gcf,'figures/exFig1_2_new.eps','Colorspace','rgb','Resolution',600);
%% extended data figure 1
Curvature = ft_read_cifti('data/templates/HCP_S1200_Atlas_Z4_pkXDZ/S1200.curvature_MSMAll.32k_fs_LR.dscalar.nii');
flatSurf = gifti('data/templates/HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.flat.32k_fs_LR.surf.gii');
midSurf = gifti('data/templates/HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii');


load('data/templates/HCP_S1200_Atlas_Z4_pkXDZ/curvature_mesh.mat');
figure;
imagesc_brainimg(curv_L_mesh,mask,0);
colormap(bone);
exportgraphics(gcf,'figures/exFig1_3.eps','Colorspace','rgb','Resolution',600);

%% 
xi = -250:2:250;
yi = -200:2:200;
xi = xi(100);
yi = yi(100);
[~,ind] = min(sum((flatSurf.vertices(:,1:2)-[xi yi]).^2,2));

figure;
brain = patch('Faces',midSurf.faces,'Vertices',midSurf.vertices);
set(brain,'facevertexcdata',Curvature.curvature_msmall(Curvature.brainstructure==1),'facecolor','interp','edgecolor','none');
colormap([0.5,0.5,0.5;bone]);
view([-6 1 1]);
ax=gca;ax.Visible='off';
axis equal;
hold on;
scatter3(midSurf.vertices(ind,1),midSurf.vertices(ind,2),midSurf.vertices(ind,3),[],'b');
exportgraphics(gcf,'figures/exFig1_4.eps','Colorspace','rgb','Resolution',600);

%% figure 2 top
analyticBrainImg(repmat(mask,1,1,T+1)==0)=-100;
t = 1;
figure('color','w');
h=imagesc_brainimg(real(analyticBrainImg(:,:,t)),mask,0,[-3 3]);
alpha(h,0.2);
ax=gca;
ax.Color=[0.97 0.97 0.97];
colormap([[0.97 0.97 0.97];jet]);
view([1 0.8 0.6]);

%%
for t = 10:11
    figure('color','w');
    h=imagesc_brainimg(real(analyticBrainImg(:,:,t)),mask,0,[-3 3]);
    alpha(h,mask);
    ax=gca;
    ax.Color=[0.97 0.97 0.97];
    colormap(jet);
    view([1 0.8 0.6]);
    figure('color','w');
    h=imagesc_brainimg(imag(analyticBrainImg(:,:,t)),mask,0,[-3 3]);
    alpha(h,mask);
    ax=gca;
    ax.Color=[0.97 0.97 0.97];
    colormap(jet);
    view([1 0.8 0.6]);
end
%% figure 2 middle
Scale = 1.5;
t = 10;
figure('color','w');
q = quiver_dsr(Scale*Ux(:,:,t),Scale*Uy(:,:,t),3,0);
q.LineWidth = 0.8;
set(q,'color','k');
axis equal; axis tight
ax=gca;ax.Visible='off';
hold on;
rectangle('Position', [139, 95, 41, 36], 'EdgeColor', 'r', 'LineWidth', 2);

%% figure 2 bottom
figure('color','w');
tiledlayout(1,6,"TileSpacing","compact","Padding","compact");
scale=1.2;
for t = 69:3:84
    nexttile
    h=imagesc_brainimg(BrainImg(:,:,t),mask,0,[MeanImg-3*StdImg MeanImg+3*StdImg]);
    alpha(h,mask);
    colormap(jet);
    hold on;
    q = quiver_dsr(scale*Ux(:,:,t),scale*Uy(:,:,t),2,0);
    q.LineWidth=1;
    q.Color='k';
    hold off;
    ax=gca;ax.XColor='none';ax.YColor='none';
    axis tight; axis equal;
    axis([139 180 95 131])
    title(['t = ' num2str(t*0.72) ' s'],'FontSize',15);
end
colorbar();
%% figure 3
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
%% figure 3a left
t = 1;

figure('color','w');
ax1=axes;
q = quiver_dsr(2*Ux(:,:,t),2*Uy(:,:,t),3,0);
q.LineWidth = 0.8;
set(q,'Color','k');
ax1.Visible='off';
axis([0 250 0 200]);axis equal;axis tight;

%% figure 3a right
figure('color','w');

ax3 = axes;
imagesc_brainimg(V_mask,1-V_mask,0,[0 1]);colormap([[1 1 1]; jet]);
ax3.Color='k';
% ax3.Visible='off';

ax1 = axes;
imagesc_brainimg(CurlMap(:,:,t),~isnan(CurlMap(:,:,t)),0,[-3 3]);
alpha(abs(CurlMap(:,:,t))/3);
colormap(ax1,winter);
ax1.Visible='off';
axis([0 250 0 200]);axis equal;axis tight;

ax2 = axes;
imagesc_brainimg(DivMap(:,:,t),~isnan(DivMap(:,:,t)),0,[-3 3]);
alpha(abs(DivMap(:,:,t))/3);
colormap(ax2,autumn);
ax2.Visible='off';
ax2.GridLineStyle='none';
axis([0 250 0 200]);axis equal;axis tight;

hold on;
q = streamslice(Ux(:,:,t),Uy(:,:,t),9);
set(q,'Color',[0.9 0.9 0.9]);
hold off

%% 3b
figure('color','w');
ha = tight_subplot(1,6,[.02 .02],[.02 .02],[.02 .02]);
Lim = [90 160 90 180];
[xi,yi]=meshgrid(Lim(1):Lim(2),Lim(3):Lim(4));
for t = 70:2:80
    axes(ha((t-68)/2));
    ax1=gca;
    axis tight;axis equal;
    ax.XTick=[];ax.YTick=[];ax.Color='none';ax.XColor='none';ax.YColor='none';ax.Color='none';
    ax1.Color='k';

    axis(Lim);
    ax3 = axes;
    imagesc_brainimg(CurlMap(:,:,t),~isnan(CurlMap(:,:,t)),0,[-3 3]);
    alpha(abs(CurlMap(:,:,t))/3);
    colormap(ax3,winter);
    axis off
    axis([0 250 0 200]);axis equal;axis tight;
    set(ax3,'Position',ax1.Position);
    axis(Lim);
    ax2 = axes;
    imagesc_brainimg(DivMap(:,:,t),~isnan(DivMap(:,:,t)),0,[-2.5 2.5]);
    alpha(abs(DivMap(:,:,t))/3);
    colormap(ax2,autumn);
    axis off
    axis([0 250 0 200]);axis equal;axis tight;
    set(ax2,'Position',ax1.Position);
    axis(Lim);

    hold on;
    q = streamslice(xi,yi,Ux(Lim(3):Lim(4),Lim(1):Lim(2),t),Uy(Lim(3):Lim(4),Lim(1):Lim(2),t),2);
    set(q,'Color',[0.9 0.9 0.9]);
    hold off
    % title(['t = ' num2str(t*0.72) ' s'],'FontSize',15);
end
set(gcf,'Position',[83,335,1277,313]);