function [DivMap, CurlMap] = meg_velocity_derivatives(Ux, Uy, V_mask, filtSigma, do_zscore, do_binarize, n_std)
%MEG_VELOCITY_DERIVATIVES Divergence/curl maps following the Project 01 convention.
%
% The optional smoothing, z-score, and binarization arguments match the
% Project 01 function. In the session-level MEG runner, z-scoring is done
% after a first pass over all chunks so thresholds are session-global.

if nargin < 5, do_zscore = 0; end
if nargin < 6, do_binarize = 0; end
if nargin < 7, n_std = 2; end
if nargin < 4 || isempty(filtSigma)
    filtSigma = 0;
end
if filtSigma < 0
    error('MEG:DivCurlSigma', 'filtSigma must be nonnegative.');
elseif filtSigma > 0
    filtWidth = round(3*filtSigma);
    imageFilter = fspecial('gaussian', filtWidth, filtSigma);
end

T = size(Ux,3);
DivMap = zeros(size(Ux));
CurlMap = zeros(size(Ux));
for t = T:-1:1
    U = sqrt(Ux(:,:,t).^2 + Uy(:,:,t).^2);
    div_oneframe = divergence(Ux(:,:,t), Uy(:,:,t));
    div_oneframe(V_mask == 0) = NaN;
    [curl_oneframe,~] = curl(Ux(:,:,t)./U, Uy(:,:,t)./U);
    curl_oneframe(V_mask == 0) = NaN;
    if filtSigma > 0
        DivMap(:,:,t) = meg_nanconv2(div_oneframe, imageFilter);
        CurlMap(:,:,t) = meg_nanconv2(curl_oneframe, imageFilter);
    else
        DivMap(:,:,t) = div_oneframe;
        CurlMap(:,:,t) = curl_oneframe;
    end
end

valid = repmat(V_mask ~= 0, 1, 1, T);
DivMap(~valid) = NaN;
CurlMap(~valid) = NaN;

if do_zscore
    divMean = mean(DivMap, 'all', 'omitnan');
    divStd = std(DivMap, [], 'all', 'omitnan');
    curlMean = mean(CurlMap, 'all', 'omitnan');
    curlStd = std(CurlMap, [], 'all', 'omitnan');
    DivMap = (DivMap - divMean) ./ divStd;
    CurlMap = (CurlMap - curlMean) ./ curlStd;
end

if do_binarize
    if do_zscore
        DivMap(DivMap < n_std & DivMap > -n_std) = 0;
        CurlMap(CurlMap < n_std & CurlMap > -n_std) = 0;
    else
        divMean = mean(DivMap, 'all', 'omitnan');
        divStd = std(DivMap, [], 'all', 'omitnan');
        curlMean = mean(CurlMap, 'all', 'omitnan');
        curlStd = std(CurlMap, [], 'all', 'omitnan');
        DivMap(DivMap < divMean+n_std*divStd & DivMap > divMean-n_std*divStd) = 0;
        CurlMap(CurlMap < curlMean+n_std*curlStd & CurlMap > curlMean-n_std*curlStd) = 0;
    end
    DivMap(DivMap > 0) = 1;
    DivMap(DivMap < 0) = -1;
    CurlMap(CurlMap > 0) = 1;
    CurlMap(CurlMap < 0) = -1;
end
end

function output = meg_nanconv2(input, kernel)
valid = isfinite(input);
values = input;
values(~valid) = 0;
if any(~valid(:))
    flat = conv2(double(valid), kernel, 'same');
    flat = flat ./ conv2(ones(size(input)), kernel, 'same');
else
    flat = ones(size(input));
end
output = conv2(values, kernel, 'same') ./ flat;
end
