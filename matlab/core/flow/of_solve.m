function [Fx,Fy] = of_solve(BrainImg,V_mask,para_model)
%OPTICALFLOW calculate the velocity vector fields (VVFs) via optical flow
%algorithm
%
% Input:
%   BrainImg: the BOLD-fMRI signal ((M+1)*(N+1)*T tensor y-, x-, 
%             t-direction)   
%   mask0: a mask for image
%   para_model: model parameters including ...
%       rho: the penalty parameter for the smooth regularization
%       tau: the penalty parameter for the uniqueness of the optimal 
%            solution
%       sigma: the parameter for detecting VVFs at certain spatial scale
%
% Output:
%   Fx,Fy: M*N*T tensors of the x-, y-components of the velocity field of 
%          original signal
%
%   by Weiyang Ding @Fudan November 4, 2023

filtSigma = para_model.sigma;
filtWidth = ceil(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

V_mask = conv2(V_mask,ones(2).*0.25,'valid');
V_mask(V_mask~=1) = NaN;
V_mask = nanconv(V_mask,imageFilter);
ind0 = find(~isnan(V_mask));
n0 = length(ind0);

V_mask0 = conv2(V_mask,ones(2).*0.25,'same');
V_mask0(1,:) = NaN;
V_mask0(end,:) = NaN;
V_mask0(:,1) = NaN;
V_mask0(:,end) = NaN;
ind1 = find(~isnan(V_mask0));
n1 = length(ind1);

[M,N,T] = size(BrainImg);% M - y-axis, N - x-axis
M = M-1; N = N-1; T = T-1;
% Gx = sparse([1:n1,1:n1]',[ind1-M;ind1+M],[-ones(n1,1).*0.5;ones(n1,1).*0.5],n1,M*N);
Gx = sparse([1:n1,1:n1]',[ind1-M;ind1],[-ones(n1,1);ones(n1,1)],n1,M*N);
Gx = Gx(:,ind0);
% Gy = sparse([1:n1,1:n1]',[ind1-1;ind1+1],[-ones(n1,1).*0.5;ones(n1,1).*0.5],n1,M*N);
Gy = sparse([1:n1,1:n1]',[ind1-1;ind1],[-ones(n1,1);ones(n1,1)],n1,M*N);
Gy = Gy(:,ind0);
GG = (Gx'*Gx + Gy'*Gy).*para_model.rho;

Fx = zeros([M,N,T]); Fy = zeros([M,N,T]);
for t = 1:T
    Dx = (BrainImg(1:M,2:N+1,t) - BrainImg(1:M,1:N,t) ...
        + BrainImg(2:M+1,2:N+1,t) - BrainImg(2:M+1,1:N,t) ...
        + BrainImg(1:M,2:N+1,t+1) - BrainImg(1:M,1:N,t+1) ...
        + BrainImg(2:M+1,2:N+1,t+1) - BrainImg(2:M+1,1:N,t+1))./4;% d/dx
    Dy = (BrainImg(2:M+1,1:N,t) - BrainImg(1:M,1:N,t) ...
        + BrainImg(2:M+1,2:N+1,t) - BrainImg(1:M,2:N+1,t) ...
        + BrainImg(2:M+1,1:N,t+1) - BrainImg(1:M,1:N,t+1) ...
        + BrainImg(2:M+1,2:N+1,t+1) - BrainImg(1:M,2:N+1,t+1))./4;% d/dy
    Dt = (BrainImg(1:M,1:N,t+1) - BrainImg(1:M,1:N,t) ...
        + BrainImg(2:M+1,1:N,t+1) - BrainImg(2:M+1,1:N,t) ...
        + BrainImg(1:M,2:N+1,t+1) - BrainImg(1:M,2:N+1,t) ...
        + BrainImg(2:M+1,2:N+1,t+1) - BrainImg(2:M+1,2:N+1,t))./4;% d/dt    
    Dx(isnan(V_mask)) = nan;Dy(isnan(V_mask)) = nan;Dt(isnan(V_mask)) = nan;
    Dx = nanconv(Dx,imageFilter);Dy = nanconv(Dy,imageFilter);Dt = nanconv(Dt,imageFilter);
    dx = Dx(ind0); dy = Dy(ind0); dt = Dt(ind0);
    b = [dx.*dt;dy.*dt];
    dt = dx.*dy;
    dx = dx.*dx+para_model.tau;
    dy = dy.*dy+para_model.tau;
    A = [spdiags(dx,0,n0,n0)+GG,spdiags(dt,0,n0,n0);...
        spdiags(dt,0,n0,n0),spdiags(dy,0,n0,n0)+GG];
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

