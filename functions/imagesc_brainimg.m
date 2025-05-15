function h = imagesc_brainimg(Map,mask,boundary,caxis_para)
Map(mask==0)=nan;
h = imagesc(Map);
set(h,'alphadata',mask==1);
ax = gca;ax.YDir = 'normal';
axis tight;axis equal;
if nargin >2 && boundary == 0
    ax.XTick=[];ax.YTick=[];ax.Color='none';ax.XColor='none';ax.YColor='none';ax.Color='none';
end
if nargin < 4
    clim([mean(Map,'all','omitnan')-3*std(Map,[],'all','omitnan') ...
        mean(Map,'all','omitnan')+3*std(Map,[],'all','omitnan')]);
else
    clim(caxis_para);
end