clear;clc
addpath(genpath('functions'))

M = 100;N = 50;T = 100;
Movie = zeros([M N T]);
%% wave
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
%% original signal
mask = ones([M N]);
[Ux_OF,Uy_OF] = velocity_field_ell2(Movie,mask,para_model);
%% hilbert
mask = ones([M N]);
Movie = reshape(Movie,[],T);
ImagPart = imag(hilbert(Movie'));
ImagPart = ImagPart';
ImagPart = reshape(ImagPart,M,N,T);
Movie = reshape(Movie,M,N,T);
BrainImg{1} = Movie;
BrainImg{2} = ImagPart;
[Ux_COF,Uy_COF] = velocity_field_multidim_ell2(BrainImg,mask,para_model);

filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
for t = T-1:-1:1
    div_oneframe = divergence(Ux_COF(:,:,t),Uy_COF(:,:,t));
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
end
%% phase gradient
mask = ones([M N]);
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
Ux_PG = -vPhaseX;
Uy_PG = -vPhaseY;
clearvars vPhaseX  vPhaseY phaseSig 

Ux_PG = convn(Ux_PG,ones(2,2,2)/8,"valid");
Uy_PG = convn(Uy_PG,ones(2,2,2)/8,"valid");
%%
[x,y] = meshgrid(1:N-1,1:M-1);
x=x+0.5;y=y+0.5;
DSR = 3;
Scale = 1;

F(80) = struct('cdata',[],'colormap',[]);
figure('color','w');
for t = 11:T-10
    clf
    tiledlayout(1,4,"TileSpacing","compact","Padding","compact");

    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    hold on
    quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux_OF(1:DSR:end,1:DSR:end,t),Scale*Uy_OF(1:DSR:end,1:DSR:end,t),0,'r');
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    title('Optical flow');

    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    hold on
    quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        50*Ux_PG(1:DSR:end,1:DSR:end,t),50*Uy_PG(1:DSR:end,1:DSR:end,t),0,'r');
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    title('Phase Gradient');

    nexttile
    imagesc(Movie(:,:,t));axis equal;axis([0.5 N+0.5 0.5 M+0.5]);
    clim([min(Movie,[],'all') max(Movie,[],'all')]);
    cc1 = colorbar();
    hold on
    quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux_COF(1:DSR:end,1:DSR:end,t),Scale*Uy_COF(1:DSR:end,1:DSR:end,t),0,'r');
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    title('Complex optical flow');

    nexttile
    imagesc(DivMap(:,:,t));axis equal;axis([1 N 1 M]);
    clim([-1 1]);
    colormap(gca,jet)
    colorbar()
    hold on
    quiver(x(1:DSR:end,1:DSR:end),y(1:DSR:end,1:DSR:end),...
        Scale*Ux_COF(1:DSR:end,1:DSR:end,t),Scale*Uy_COF(1:DSR:end,1:DSR:end,t),0,'r');
    hold off
    ax=gca;ax.XColor='none';ax.YColor='none';axis equal;axis tight;
    title('Complex optical flow');

    sgtitle(['t = ' num2str(t-10)],"FontSize",15);
    set(gcf,'Position',[234.5,134.5,944,384]);

    drawnow;
    F(t-10) = getframe(gcf);
end

createmovie(F,'data/videos/SupplemantaryVideo1',3);