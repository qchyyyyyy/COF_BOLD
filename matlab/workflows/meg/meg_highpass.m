function [filtered, X0] = meg_highpass(data, TR, cutoffSeconds)
%MEG_HIGHPASS SPM-compatible DCT high-pass residual forming matrix.
%
% Equivalent to spm_filter for one contiguous session with K.row=1:T.

if ndims(data) ~= 2 || TR <= 0 || cutoffSeconds <= 0
    error('MEG:WorkflowHighpassArguments', 'Expected time-by-variable data and positive timing.');
end
nTime = size(data,1);
nDct = fix(2 * nTime * TR / cutoffSeconds + 1);
if nDct <= 1
    X0 = zeros(nTime,0);
    filtered = data;
    return
end
X = dct_basis(nTime, nDct);
X0 = X(:,2:end);
filtered = data - X0 * (X0' * data);
end

function X = dct_basis(n, k)
% Columns match spm_dctmtx(N,K): constant followed by orthonormal cosines.
t = (0:n-1)' + 0.5;
X = zeros(n,k);
X(:,1) = 1 ./ sqrt(n);
for j = 2:k
    X(:,j) = sqrt(2/n) .* cos(pi .* (j-1) .* t ./ n);
end
end
