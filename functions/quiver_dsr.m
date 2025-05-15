function qq = quiver_dsr(Ux,Uy,dsr,quiver_para)
V_mask = zeros(size(Ux));
V_mask(1:dsr:end,1:dsr:end) = 1;
V_mask(Ux==0) = 0;
[M,N] = size(Ux);
[x,y] = meshgrid(1:N,1:M);
if isempty(dsr) || nargin < 3
    dsr = 1;
end
if nargin < 4
    quiver_para = 0;
end
qq = quiver(x(V_mask==1),y(V_mask==1),...
    Ux(V_mask==1),Uy(V_mask==1),...
    quiver_para);
end
