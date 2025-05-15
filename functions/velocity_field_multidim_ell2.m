function [Fx,Fy] = velocity_field_multidim_ell2(BrainImg,mask,para_model)
% Input:
%   BrainImg: an 1*2 cell that contains the real and imaginary parts of 
%       the analytic BOLD-fMRI signal ((M+1)*(N+1)*T tensor y-, x-, t-direction)   
%   mask0: a mask for image
%   para_model: model parameters including ...
%       rho: the penalty parameter for the smooth regularization
%       tau: the penalty parameter for the uniqueness of the optimal solution
%       sigma: the parameter for detecting VVFs at certain spatial scale
%
% Output:
%   Fx,Fy: M*N*T tensors of the x-, y-components of the velocity field of analytic signal
%   by Weiyang Ding @Fudan November 4, 2023

filtSigma = para_model.sigma;
filtWidth = ceil(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

V_mask = conv2(mask,ones(2).*0.25,'valid');
V_mask(V_mask==0) = NaN;
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

[M,N,T] = size(BrainImg{1});% M - y-axis, N - x-axis
M = M-1; N = N-1; T = T-1;
Gx = sparse([1:n1,1:n1]',[ind1-M;ind1],[-ones(n1,1);ones(n1,1)],n1,M*N);
Gx = Gx(:,ind0);
Gy = sparse([1:n1,1:n1]',[ind1-1;ind1],[-ones(n1,1);ones(n1,1)],n1,M*N);
Gy = Gy(:,ind0);
GG = (Gx'*Gx + Gy'*Gy).*para_model.rho;

Fx = zeros([M,N,T]); Fy = zeros([M,N,T]);
for t = 1:T
    b = zeros(n0*2,1);
    wt = zeros(n0,1);
    wx = zeros(n0,1);
    wy = zeros(n0,1);
    for k = 1:K
        Dx = (BrainImg{k}(1:M,2:N+1,t) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t) - BrainImg{k}(2:M+1,1:N,t) ...
            + BrainImg{k}(1:M,2:N+1,t+1) - BrainImg{k}(1:M,1:N,t+1) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(2:M+1,1:N,t+1))./4;% d/dx
        Dy = (BrainImg{k}(2:M+1,1:N,t) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t) - BrainImg{k}(1:M,2:N+1,t) ...
            + BrainImg{k}(2:M+1,1:N,t+1) - BrainImg{k}(1:M,1:N,t+1) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(1:M,2:N+1,t+1))./4;% d/dy
        Dt = (BrainImg{k}(1:M,1:N,t+1) - BrainImg{k}(1:M,1:N,t) ...
            + BrainImg{k}(2:M+1,1:N,t+1) - BrainImg{k}(2:M+1,1:N,t) ...
            + BrainImg{k}(1:M,2:N+1,t+1) - BrainImg{k}(1:M,2:N+1,t) ...
            + BrainImg{k}(2:M+1,2:N+1,t+1) - BrainImg{k}(2:M+1,2:N+1,t))./4;% d/dt
        Dx(isnan(V_mask)) = nan;Dy(isnan(V_mask)) = nan;Dt(isnan(V_mask)) = nan;
        Dx = nanconv(Dx,imageFilter);Dy = nanconv(Dy,imageFilter);Dt = nanconv(Dt,imageFilter);
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
    A = [spdiags(wx,0,n0,n0)+GG,spdiags(wt,0,n0,n0);...
        spdiags(wt,0,n0,n0),spdiags(wy,0,n0,n0)+GG];
    u = -A\b;
    Ux = zeros(M,N);
    Uy = zeros(M,N);
    Ux(ind0) = u(1:n0);
    Uy(ind0) = u(n0+1:n0*2);
    Fx(:,:,t) = Ux;
    Fy(:,:,t) = Uy;
    %disp([num2str(t),'-th step finishes'])
end

end

