function [Fx,Fy] = meg_cof_solve(BrainImg,mask,para_model)
%MEG_COF_SOLVE Project 01 COF solver with explicit sigma=0 support.
%
% BrainImg is a cell array containing the real and imaginary parts of the
% analytic signal. sigma=0 means no Gaussian smoothing (identity kernel).

% Based on the Project 01 implementation by Weiyang Ding, 2023.

filtSigma = para_model.sigma;
if filtSigma < 0
    error('MEG:COFSigma', 'sigma must be nonnegative.');
elseif filtSigma == 0
    imageFilter = 1;
else
    filtWidth = ceil(3*filtSigma);
    imageFilter = fspecial('gaussian',filtWidth,filtSigma);
end

V_mask = conv2(mask,ones(2).*0.25,'valid');
V_mask(V_mask~=1) = NaN;
if filtSigma > 0
    V_mask = nanconv(V_mask,imageFilter);
end
ind0 = find(~isnan(V_mask));
n0 = length(ind0);

V_mask0 = conv2(V_mask,ones(2).*0.25,'same');
V_mask0(1,:) = NaN;
V_mask0(end,:) = NaN;
V_mask0(:,1) = NaN;
V_mask0(:,end) = NaN;
ind1 = find(~isnan(V_mask0));
n1 = length(ind1);

K = length(BrainImg);
[M,N,T] = size(BrainImg{1});
M = M-1; N = N-1; T = T-1;
Gx = sparse([1:n1,1:n1]',[ind1-M;ind1], ...
    [-ones(n1,1);ones(n1,1)],n1,M*N);
Gx = Gx(:,ind0);
Gy = sparse([1:n1,1:n1]',[ind1-1;ind1], ...
    [-ones(n1,1);ones(n1,1)],n1,M*N);
Gy = Gy(:,ind0);
GG = (Gx'*Gx + Gy'*Gy).*para_model.rho;

Fx = zeros([M,N,T]);
Fy = zeros([M,N,T]);
V_mask = conv2(mask,ones(2).*0.25,'valid');
V_mask(V_mask~=1) = NaN;
for t = 1:T
    b = zeros(n0*2,1);
    wt = zeros(n0,1);
    wx = zeros(n0,1);
    wy = zeros(n0,1);
    for k = 1:K
        Dx = (BrainImg{k}(1:M,2:N+1,t) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t) - BrainImg{k}(2:M+1,1:N,t) ...
            + BrainImg{k}(1:M,2:N+1,t+1) - BrainImg{k}(1:M,1:N,t+1) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(2:M+1,1:N,t+1))./4;
        Dy = (BrainImg{k}(2:M+1,1:N,t) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t) - BrainImg{k}(1:M,2:N+1,t) ...
            + BrainImg{k}(2:M+1,1:N,t+1) - BrainImg{k}(1:M,1:N,t+1) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(1:M,2:N+1,t+1))./4;
        Dt = (BrainImg{k}(1:M,1:N,t+1) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,1:N,t+1) - BrainImg{k}(2:M+1,1:N,t) ...
            + BrainImg{k}(1:M,2:N+1,t+1) - BrainImg{k}(1:M,2:N+1,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(2:M+1,2:N+1,t))./4;
        Dx(isnan(V_mask)) = NaN;
        Dy(isnan(V_mask)) = NaN;
        Dt(isnan(V_mask)) = NaN;
        if filtSigma > 0
            Dx = nanconv(Dx,imageFilter);
            Dy = nanconv(Dy,imageFilter);
            Dt = nanconv(Dt,imageFilter);
        end
        dx = Dx(ind0); dy = Dy(ind0); dt = Dt(ind0);
        b = b + [dx.*dt;dy.*dt];
        wt = wt + dx.*dy;
        wx = wx + dx.*dx;
        wy = wy + dy.*dy;
    end
    b = b./K;
    wt = wt./K;
    wx = wx./K + para_model.tau;
    wy = wy./K + para_model.tau;
    A = [spdiags(wx,0,n0,n0)+GG,spdiags(wt,0,n0,n0); ...
        spdiags(wt,0,n0,n0),spdiags(wy,0,n0,n0)+GG];
    u = -A\b;
    Ux = zeros(M,N);
    Uy = zeros(M,N);
    Ux(ind0) = u(1:n0);
    Uy(ind0) = u(n0+1:n0*2);
    Fx(:,:,t) = Ux;
    Fy(:,:,t) = Uy;
end
end
