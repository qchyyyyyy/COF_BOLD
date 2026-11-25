function [cope,varcope,tstat,beta_gls] = fit_gls_ar1(Y, X, c)
% Y: T×V  (time × vertices)
% X: T×P  (design, NO HRF convolution; include intercept)
% c: P×1  (contrast, e.g., [0; +1; -1] for 2-back > 0-back)

[~,V] = size(Y);
P = size(X,2);

% ---- Step 1: OLS per vertex ----
XTX = X' * X;
XTY = X' * Y;
beta_ols = XTX \ XTY;              % P×V
E = Y - X * beta_ols;              % residuals T×V

% ---- Step 2: AR(1) rho per vertex (lag-1 autocorr / Yule–Walker) ----
E1 = E(1:end-1,:); E2 = E(2:end,:);
num = sum(E1 .* E2, 1);
den = sum(E1 .* E1, 1) + eps;
rho_hat = num ./ den;
rho_hat = max(min(rho_hat, 0.99), -0.99);   % stabilize |rho|<1

% ---- Step 3: Prais–Winsten whitening & GLS (OLS on whitened) ----
cope     = zeros(1,V);
varcope  = zeros(1,V);
tstat    = zeros(1,V);
beta_gls = zeros(P,V);

for v = 1:V
    rho = rho_hat(v);

    % Prais–Winsten transform (keeps the first observation)
    % y_w(1) = sqrt(1-rho^2)*y(1)
    % y_w(t) = y(t) - rho*y(t-1), t=2..T
    yw = [sqrt(1-rho^2)*Y(1,v);  Y(2:end,v) - rho*Y(1:end-1,v)];
    Xw = [sqrt(1-rho^2)*X(1,:);  X(2:end,:) - rho*X(1:end-1,:)];

    % GLS == OLS on whitened data
    XtX = Xw' * Xw;
    XtY = Xw' * yw;
    b   = XtX \ XtY;                 % P×1
    beta_gls(:,v) = b;

    yhat = Xw * b;
    resw = yw - yhat;
    dof  = size(Xw,1) - rank(Xw);
    sig2 = (resw' * resw) / max(dof,1);      % scalar

    c_cov = (c' / XtX) * c;                  % scalar
    cope(1,v)    = (c' * b);
    varcope(1,v) = sig2 * c_cov;
    tstat(1,v)   = cope(1,v) / sqrt(varcope(1,v) + eps);
end
end