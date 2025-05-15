function [Ux,Uy] = recon_velfield(DivMap,CurlMap,V_mask,lambda,rho)
[M,N] = size(V_mask);
ind0 = find(V_mask==1);
n0 = length(ind0);
xs = zeros(size(ind0));xa=xs;xsv=xs;xav=xs;ys=xs;ya=xs;ysv=xs;yav=xs;
for node = n0:-1:1
    [y,x] = ind2sub(size(V_mask),ind0(node));
    if V_mask(y,x+1)==1 && V_mask(y,x-1)==1
        xs(node)=-M;xa(node)=M;
        xsv(node) = -1/2;xav(node) = 1/2;
    elseif V_mask(y,x+1)==0
        xs(node)=-M;xa(node)=0;
        xsv(node) = -1;xav(node) = 1;
    elseif V_mask(y,x-1)==0
        xs(node)=0;xa(node)=M;
        xsv(node) = -1;xav(node) = 1;
    end
    if V_mask(y+1,x)==1 && V_mask(y-1,x)==1
        ys(node)=-1;ya(node)=1;
        ysv(node) = -1/2;yav(node) = 1/2;
    elseif V_mask(y+1,x)==0
        ys(node)=-1;ya(node)=0;
        ysv(node) = -1;yav(node) = 1;
    elseif V_mask(y-1,x)==0
        ys(node)=0;ya(node)=1;
        ysv(node) = -1;yav(node) = 1;
    end
end
Gx = sparse([1:n0,1:n0]',[ind0+xs;ind0+xa],[xsv;xav],n0,M*N);
Gx = Gx(:,ind0);
Gy = sparse([1:n0,1:n0]',[ind0+ys;ind0+ya],[ysv;yav],n0,M*N);
Gy = Gy(:,ind0);

% Gx = sparse([1:n0,1:n0]',[ind0-M;ind0+M],[-ones(n0,1)/2;ones(n0,1)/2],n0,M*N);
% Gx = Gx(:,ind0);
% Gy = sparse([1:n0,1:n0]',[ind0-1;ind0+1],[-ones(n0,1)/2;ones(n0,1)/2],n0,M*N);
% Gy = Gy(:,ind0);
if nargin < 4
    lambda = 1e-6;
    rho = 1+1e-2;
end
c = CurlMap(ind0);
d = DivMap(ind0);
GTG = [(Gx'*Gx+Gy'*Gy).*rho,Gx'*Gy-Gy'*Gx;Gy'*Gx-Gx'*Gy,(Gx'*Gx+Gy'*Gy).*rho];
GTb = [-Gy'*c+Gx'*d;Gx'*c+Gy'*d];

u = (GTG+lambda*speye(2*n0))\GTb;

Ux = vec2map(u(1:end/2),V_mask);
Uy = vec2map(u((end/2)+1:end),V_mask);