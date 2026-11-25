function hrf = meg_hrf(RT, p, microtimeResolution)
%MEG_HRF Dependency-light numerical equivalent of SPM spm_hrf.

if nargin < 2 || isempty(p), p = [6 16 1 1 6 0 32]; end
if nargin < 3 || isempty(microtimeResolution), microtimeResolution = 16; end
if RT <= 0 || numel(p) ~= 7 || microtimeResolution < 1
    error('MEG:WorkflowHRFArguments', 'Invalid HRF arguments.');
end
p = double(p(:)');
dt = RT / microtimeResolution;
u = (0:ceil(p(7)/dt)) - p(6)/dt;
hrf = gamma_rate_pdf(u, p(1)/p(3), dt/p(3)) - ...
    gamma_rate_pdf(u, p(2)/p(4), dt/p(4)) ./ p(5);
indices = (0:floor(p(7)/RT)) .* microtimeResolution + 1;
hrf = hrf(indices)';
normalizer = sum(hrf);
if ~isfinite(normalizer) || abs(normalizer) <= eps
    error('MEG:WorkflowHRFNormalization', 'HRF normalization is invalid.');
end
hrf = hrf ./ normalizer;
end

function y = gamma_rate_pdf(x, shape, rate)
y = zeros(size(x));
positive = x > 0;
xp = x(positive);
y(positive) = exp((shape-1).*log(xp) + shape.*log(rate) - ...
    rate.*xp - gammaln(shape));
if shape == 1, y(x == 0) = rate; end
end
