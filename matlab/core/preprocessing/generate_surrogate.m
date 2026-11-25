function is_valid = generate_surrogate(input_file, output_dir, surrogateMethod, cii_mask_left, cii_mask_right, randomSeed)
%GENERATE_SURROGATE Generate phase-randomized surrogate data from preprocessed 2-D fMRI.
%
% Input:
%   input_file       Preprocessed TC_mesh_r_f_zscore.mat file.
%   output_dir       Directory for the generated surrogate file.
%   surrogateMethod  '2D' or '3D'. Default is '2D'.
%
% Output:
%   output_dir/TC_surrogate.mat with data.BrainImg_lh and data.BrainImg_rh.

if nargin < 3 || isempty(surrogateMethod)
    surrogateMethod = '2D';
end

if nargin >= 6 && ~isempty(randomSeed)
    rng(randomSeed);
end

is_valid = false;

if ~exist(input_file, 'file')
    disp([input_file ' not exist']);
    return
end

if nargin < 5 || isempty(cii_mask_left) || isempty(cii_mask_right)
    load(template_path('coordinate_cii_v2.mat'), 'cii_mask_left', 'cii_mask_right');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_file = fullfile(output_dir, 'TC_surrogate.mat');
if exist(output_file, 'file')
    disp([output_file ' already exist']);
    is_valid = true;
    return
end

S = load(input_file, 'data');
data = S.data;

if ~isfield(data, 'BrainImg_lh') || ~isfield(data, 'BrainImg_rh')
    disp([input_file ' does not contain data.BrainImg_lh and data.BrainImg_rh']);
    return
end

TC_lh_mesh = data.BrainImg_lh;
TC_rh_mesh = data.BrainImg_rh;

if ~isequal(size(TC_lh_mesh, 1), size(cii_mask_left, 1)) || ~isequal(size(TC_lh_mesh, 2), size(cii_mask_left, 2))
    disp('Left hemisphere data and mask sizes do not match');
    return
end

if ~isequal(size(TC_rh_mesh, 1), size(cii_mask_right, 1)) || ~isequal(size(TC_rh_mesh, 2), size(cii_mask_right, 2))
    disp('Right hemisphere data and mask sizes do not match');
    return
end

n_left_cols = size(TC_lh_mesh, 2);
TC = cat(2, TC_lh_mesh, TC_rh_mesh);
mask = [cii_mask_left cii_mask_right];
clearvars TC_lh_mesh TC_rh_mesh

% Keep the same midline-column handling as preprocess_surrogate.m.
TC(:, n_left_cols, :) = [];
mask(:, n_left_cols) = [];

switch upper(surrogateMethod)
    case '3D'
        TC = fft_3d_surrogate(TC);
    case '2D'
        TC = fft_2d_surrogate(grid_to_vector(TC, mask));
        TC = vector_to_grid(TC, mask);
    otherwise
        error('Unknown surrogateMethod: %s. Use ''2D'' or ''3D''.', surrogateMethod);
end

data.BrainImg_lh = TC(:, 1:n_left_cols, :);
data.BrainImg_rh = TC(:, n_left_cols:end, :);
data.SurrogateMethod = surrogateMethod;
data.SurrogateSourceFile = input_file;
data.SurrogateOutputFile = output_file;

save(output_file, 'data', '-v7.3');
disp([output_file ' is successfully generated']);
is_valid = true;

end

function dataOut = fft_3d_surrogate(dataIn)

if mod(size(dataIn,3),2)==0
    dataIn = dataIn(:,:,1:end-1);
end
midPoint = floor(size(dataIn) / 2) + 1;
indMid = sub2ind(size(dataIn), midPoint(1), midPoint(2), midPoint(3));
freqData = fftshift(fftn(dataIn));
phaseData = angle(freqData);
absData = abs(freqData);
phaseDataRand = phaseData;
phaseDataRand(1:indMid-1) = rand(indMid-1, 1) .* 2 .* pi - pi;
for iInd = 1:indMid-2
    phaseDataRand(indMid+iInd) = -phaseDataRand(indMid-iInd);
end
freqDataNew = ifftshift(absData .* exp(1i .* phaseData + 1i .* phaseDataRand));
dataOut = real(ifftn(freqDataNew));

end

function dataOut = fft_2d_surrogate(dataIn)

midPoint = floor(size(dataIn, 2) / 2) + 1;
freqData = fftshift(fft(dataIn, [], 2));
phaseData = angle(freqData);
absData = abs(freqData);
phaseDataRand = phaseData;
phaseDataRand(:, 1:midPoint-1) = rand(size(phaseData, 1), midPoint-1) .* 2 .* pi - pi;
for iInd = 1:midPoint-2
    phaseDataRand(:, midPoint+iInd) = -phaseDataRand(:, midPoint-iInd);
end
freqDataNew = ifftshift(absData .* exp(1i .* phaseData + 1i .* phaseDataRand));
dataOut = real(ifft(freqDataNew, [], 2));

end
