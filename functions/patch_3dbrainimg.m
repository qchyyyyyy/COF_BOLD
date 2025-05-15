function patch_3dbrainimg(Map,templatename_3d,templatename_2d,Clims,viewdirection)
surf_3d = gifti(templatename_3d);
surf_2d = gifti(templatename_2d);
x = surf_2d.vertices(:,1);y = surf_2d.vertices(:,2);

Vortex = zeros(size(x));
for node = 1:length(x)
    if x(node)==0
        Vortex(node) = nan;
    else
        if size(Map,1)==201
            Vortex(node) = Map(round((y(node)+200)/2+1),round((x(node)+250)/2+1));
        elseif size(Map,1)==200
            Vortex(node) = Map(round((y(node)+200)/2+1.5),round((x(node)+250)/2+1.5));
        end
    end
end

% Vortex(isnan(Vortex))=0;
if nargin <4
    Clims = [min(Vortex), max(Vortex)];
end

if Clims(1)<0
    Vortex(isnan(Vortex)) = Clims(1)*1.1;
elseif Clims(1)==0
    Vortex(isnan(Vortex)) = -0.1*(Clims(2)-Clims(1));
else
    Vortex(isnan(Vortex)) = Clims(1)*0.9;
end



brain = patch('Faces',surf_3d.faces,'Vertices',surf_3d.vertices);
set(brain,'facevertexcdata',Vortex,'facecolor','interp','edgecolor','none');
colormap([0.5,0.5,0.5;jet]);

% camlight('headlight')
% material dull
% brighten(.2);
clim(Clims);
% caxis([mean(Vortex,'all','omitnan')-3*std(Vortex,[],'all','omitnan') ...
%     mean(Vortex,'all','omitnan')+3*std(Vortex,[],'all','omitnan')]);
%     colorbar;
ax=gca;
ax.XTick=[];ax.YTick=[];ax.ZTick=[];
ax.XColor='none';ax.YColor='none';ax.ZColor='none';ax.Color='none';
if nargin <5
    view([-6 1 1]);
else
    view(viewdirection);
end
axis equal;grid off;
hold off;