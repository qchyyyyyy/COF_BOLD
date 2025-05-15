clear;clc
addpath(genpath('functions'))

M = 100;N = 50;T = 100;
Movie = zeros([M N T]);
%% wave ref: https://github.com/tsb46/BOLD_WAVES/blob/master/simulation.ipynb
y=-1;
w=0.5;
t_vec = linspace(1,2*pi*20,400);
t_decay_vec = 0.5*cos(linspace(1,2*pi*2,400))+0.5;
[X, Y] = meshgrid(linspace(-6,6,N),linspace(-6,6,M));
c = exp(-(X.^2+(Y-1.5).^2)./8);
d = exp(-(X.^2+(Y+1.5).^2)./8);

for t = 1:T
    Movie(:,:,t) = 2*exp(y*t_decay_vec(t))*(cos(w*t_vec(t))*c-sin(w*t_vec(t))*d)+1e-5*randn([M N]);
end
para_model.tau = 1e-4;
para_model.rho = 1e-3;
para_model.sigma = 1e-3;

LineWidth = 0.8;
Color = 'k';

%% figure 1a (optical flow)
mask = ones([M N]);
[Ux,Uy] = velocity_field_ell2(Movie,mask,para_model);

[x,y] = meshgrid(1:N-1,1:M-1);
x=x+0.5;y=y+0.5;
DSR = 3;
Scale = 1;

figure('color','w');set(gcf,"Position",[1,10,2100,700])
set(gcf,'Position',[18,117,1402,560]);
tiledlayout(4,5,"TileSpacing","compact","Padding","compact");
for t = 22:4:38
    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    clim([-3 3]);
    % colorbar()
    hold on;
    qq = quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux(1:DSR:end,1:DSR:end,t),Scale*Uy(1:DSR:end,1:DSR:end,t),0,'r');
    qq.LineWidth = LineWidth;
    qq.Color = Color;
    hold off;
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    title(['t = ' num2str(t/2)],"FontSize",14);
end
%% figure 1b (phase gradient)
Movie = reshape(Movie,[],T);
phaseSig = angle(hilbert(Movie'));
phaseSig = phaseSig';
phaseSig = reshape(phaseSig,M,N,T);
Movie = reshape(Movie,M,N,T);
vPhaseX = zeros(size(phaseSig)) ;
vPhaseY = zeros(size(phaseSig)) ;
for iTime = 1:size(phaseSig,3)
    for iX = 1:size(phaseSig,1)
        vPhaseX(iX,2:end-1,iTime) = (anglesubtract(phaseSig(iX,3:end,iTime),phaseSig(iX,1:end-2,iTime)))/2 ;
    end
    for iY = 1:size(phaseSig,2)
        vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime),phaseSig(1:end-2,iY,iTime)))/2 ;
    end
end
Ux = -vPhaseX;
Uy = -vPhaseY;
clearvars vPhaseX  vPhaseY phaseSig 

[x,y] = meshgrid(1:N,1:M);
DSR = 3;
Scale = 70;

for t = 22:4:38
    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    clim([-3 3]);
    hold on
    qq = quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux(1:DSR:end,1:DSR:end,t),Scale*Uy(1:DSR:end,1:DSR:end,t),0);
    qq.LineWidth = LineWidth;
    qq.Color = Color;
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    % title(['t = ' num2str(t/2)],"FontSize",20);
end
%% figure 1c & 1d (complex optical flow)
mask = ones([M N]);
Movie = reshape(Movie,[],T);
ImagPart = imag(hilbert(Movie.'));
ImagPart = ImagPart.';
ImagPart = reshape(ImagPart,M,N,T);
Movie = reshape(Movie,M,N,T);
BrainImg{1} = Movie;
BrainImg{2} = ImagPart;
[Ux,Uy] = velocity_field_multidim_ell2(BrainImg,mask,para_model);

filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
for t = T-1:-1:1
    div_oneframe = divergence(Ux(:,:,t),Uy(:,:,t));
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
end

%% 1c
[x,y] = meshgrid(1:N-1,1:M-1);
x=x+0.5;y=y+0.5;
DSR = 3;
Scale = 1;

for t = 22:4:38
    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    clim([-3 3]);
    hold on
    qq = quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux(1:DSR:end,1:DSR:end,t),Scale*Uy(1:DSR:end,1:DSR:end,t),0,'r');
    qq.LineWidth = LineWidth;
    qq.Color = Color;
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    % title(['t = ' num2str(t/2)],"FontSize",20);
end
%% 1d 
for t = 22:4:38
    nexttile
    imagesc(DivMap(:,:,t));axis equal;axis([1 N 1 M]);
    clim([-1 1]);
    hold on
    qq = quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux(1:DSR:end,1:DSR:end,t),Scale*Uy(1:DSR:end,1:DSR:end,t),0,'r');
    qq.LineWidth = LineWidth;
    qq.Color = Color;
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;colormap(ax,jet)
    % title(['t = ' num2str(t/2)],"FontSize",20);
end

% colormap(gcf,'bluewhitered')