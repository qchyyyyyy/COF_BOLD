function Outputs = rest_principal_modes(cfg)
%FIG1F_I_PRINCIPAL_PATTERN_MODE Compute numerical summaries of the principal propagation mode.
%

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

cfg = rest_config(cfg);

geometry = make_mode_geometry(cfg);
displayGeometry = make_display_geometry(cfg);

Outputs = struct();
Outputs.cfg = cfg;
Outputs.geometry = geometry;

if ~exist(cfg.groupMetricDir, 'dir')
    mkdir(cfg.groupMetricDir);
end
if cfg.saveData && ~exist(cfg.figureDataDir, 'dir')
    mkdir(cfg.figureDataDir);
end

Outputs.modeFiles = compute_or_load_modes(cfg, geometry);

if cfg.makeFigures
    gradientMesh=load_gradient_mesh(cfg,geometry);
    P=struct('method',{},'timeMap',{},'gradientBin',{},'Div',{},'Curl',{},'Ux',{},'Uy',{},'div_r',{},'div_p',{},'velocity_r',{},'velocity_p',{});
    for k=1:numel(cfg.plotMethods)
        method=cfg.plotMethods{k};M=load_required_mode(cfg,method);
        [D,C]=mode_to_maps(M.u1,geometry);
        [D,C]=zscore_mode_maps(D,C,geometry.V_mask);
        [Db,Gb]=bin_by_gradient(D(geometry.V_mask==1),gradientMesh(geometry.V_mask==1),100);
        frames=valid_time_window(cfg.modeTimeWindow,numel(M.v1));
        timeMap=Db*M.v1(frames).';
        threshold=threshold_signed_map(D,geometry.V_mask,2);
        [dr,dp]=corr(threshold(geometry.V_mask==1),gradientMesh(geometry.V_mask==1),'Rows','complete');
        if strcmp(method,'PG'),D=-D;C=-C;end
        [Ux,Uy]=reconstruct_velocity(D,C,geometry.V_mask);
        [gx,gy]=gradient(gradientMesh);
        [vr,vp]=corr([Ux(geometry.V_mask==1);Uy(geometry.V_mask==1)], ...
                    [gx(geometry.V_mask==1);gy(geometry.V_mask==1)],'Rows','complete');
        P(k)=struct('method',method,'timeMap',timeMap,'gradientBin',Gb,'Div',D,'Curl',C, ...
            'Ux',Ux,'Uy',Uy,'div_r',dr,'div_p',dp,'velocity_r',vr,'velocity_p',vp);
    end
    Outputs.renderFile=write_plot_data('fig1f_i_principal_pattern_mode_render_saved', ...
        {P,geometry,displayGeometry,gradientMesh,cfg,fullfile(cfg.figureDir,'principal_modes')},6);
    Outputs.displayStatistics=P;
end


end

function geometry = make_mode_geometry(cfg)

[~, ~, ~, ~, V_mask_L] = load_cortical_atlas(cfg.modeParcellation, cfg.hemisphere(1), cfg.dataRoot);
[~, ~, ~, ~, V_mask_R] = load_cortical_atlas(cfg.modeParcellation, cfg.hemisphere(2), cfg.dataRoot);

V_mask = [V_mask_L V_mask_R];
V_mask_dsr = zeros(size(V_mask));
V_mask_dsr(1:cfg.downSampRate:end, 1:cfg.downSampRate:end) = 1;
V_mask_dsr(V_mask == 0) = 0;

[y, x] = ind2sub(size(V_mask_dsr), find(V_mask_dsr == 1));

geometry = struct();
geometry.V_mask = V_mask;
geometry.V_mask_dsr = V_mask_dsr;
geometry.x = x;
geometry.y = y;
geometry.xTarget = 1:size(V_mask, 2);
geometry.yTarget = 1:size(V_mask, 1);

end

function displayGeometry = make_display_geometry(cfg)

[mask_L, ~, oi_x_L, oi_y_L, ~] = load_cortical_atlas(cfg.parcellation, cfg.hemisphere(1), cfg.dataRoot);
[~, ~, oi_x_R, oi_y_R, ~] = load_cortical_atlas(cfg.parcellation, cfg.hemisphere(2), cfg.dataRoot);

displayGeometry = struct();
displayGeometry.oi_x = [oi_x_L; oi_x_R + size(mask_L, 2)];
displayGeometry.oi_y = [oi_y_L; oi_y_R];

end

function modeFiles = compute_or_load_modes(cfg, geometry)

modeFiles = struct();

for i_method = 1:numel(cfg.modeMethods)
    method = cfg.modeMethods{i_method};
    outFile = mode_file(cfg, method);
    modeFiles.(method) = outFile;

    if exist(outFile, 'file') && ~cfg.forceModeComputation
        fprintf('Mode file already exists for %s: %s\n', method, outFile);
        continue
    end

    subjectList = list_subjects_with_vvfs(cfg.vvfDataDir, method);
    if isempty(subjectList)
        error('No subjects with %s velocity fields found under: %s', method, cfg.vvfDataDir);
    end

    fprintf('Computing principal %s mode from %d subjects.\n', method, numel(subjectList));
    [u1, s1, v1, vi_length] = compute_one_method_mode(cfg, geometry, method, subjectList);
    save(outFile, 'u1', 's1', 'v1', 'vi_length', 'subjectList', 'method', 'cfg', '-v7.3');
end

end

function [u1, s1, v1, vi_length] = compute_one_method_mode(cfg, geometry, method, subjectList)

nPixels = nnz(geometry.V_mask_dsr);
AAT = zeros(2 * nPixels, 2 * nPixels);
vi_length = zeros(numel(subjectList), 1);

for i_sub = 1:numel(subjectList)
    subjectID = subjectList{i_sub};
    fprintf('AAT pass | %s | subject %d/%d: %s\n', method, i_sub, numel(subjectList), subjectID);
    A = subject_divcurl_matrix(cfg, geometry, subjectID, method);
    vi_length(i_sub) = size(A, 2);
    AAT = AAT + A * A.';
end

opts = struct();
opts.isreal = true;
opts.issym = true;
[u1, lambda] = eigs(AAT, 1, 'largestabs', opts);
s1 = sqrt(lambda(1, 1));
clear AAT;

scoreCells = cell(numel(subjectList), 1);
for i_sub = 1:numel(subjectList)
    subjectID = subjectList{i_sub};
    fprintf('Score pass | %s | subject %d/%d: %s\n', method, i_sub, numel(subjectList), subjectID);
    A = subject_divcurl_matrix(cfg, geometry, subjectID, method);
    scoreCells{i_sub} = (A.' * u1) ./ s1;
end

v1 = vertcat(scoreCells{:});
vi_length = cellfun(@numel, scoreCells);

end

function A = subject_divcurl_matrix(cfg, geometry, subjectID, method)

subjectDir = fullfile(cfg.vvfDataDir, subjectID);
[DivMap, CurlMap] = load_subject_divcurl(cfg, geometry.V_mask, subjectDir, method);

DivMap = grid_to_vector(DivMap, geometry.V_mask_dsr);
CurlMap = grid_to_vector(CurlMap, geometry.V_mask_dsr);
A = [DivMap; CurlMap];

end

function [DivMap, CurlMap] = load_subject_divcurl(cfg, V_mask, subjectDir, method)

leftFile = fullfile(subjectDir, ['VelField_lh_' method '.mat']);
rightFile = fullfile(subjectDir, ['VelField_rh_' method '.mat']);
if ~exist(leftFile, 'file') || ~exist(rightFile, 'file')
    error('Missing velocity-field files for %s under %s.', method, subjectDir);
end

L = load(leftFile, 'VelField');
R = load(rightFile, 'VelField');

if strcmp(method, 'PG')
    Ux = cat(2, L.VelField.Ux(:, 1:end-1, :), R.VelField.Ux);
    Uy = cat(2, L.VelField.Uy(:, 1:end-1, :), R.VelField.Uy);
    Ux = convn(Ux, ones(2, 2, 2) ./ 8, 'valid');
    Uy = convn(Uy, ones(2, 2, 2) ./ 8, 'valid');
else
    Ux = cat(2, L.VelField.Ux, R.VelField.Ux);
    Uy = cat(2, L.VelField.Uy, R.VelField.Uy);
end

[DivMap, CurlMap] = velocity_derivatives(Ux, Uy, V_mask, cfg.filtSigma, cfg.do_zscore, cfg.do_binarize);

end

function subjectList = list_subjects_with_vvfs(dataDir, method)

if ~exist(dataDir, 'dir')
    subjectList = {};
    return
end

items = dir(dataDir);
items = items([items.isdir]);
subjectList = {};
for i_item = 1:numel(items)
    name = items(i_item).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue
    end
    if exist(fullfile(dataDir, name, ['VelField_lh_' method '.mat']), 'file') && ...
            exist(fullfile(dataDir, name, ['VelField_rh_' method '.mat']), 'file')
        subjectList{end+1, 1} = name; %#ok<AGROW>
    end
end

subjectList = sort(subjectList);

end

function outFile = mode_file(cfg, method)

outFile = fullfile(cfg.groupMetricDir, ...
    [cfg.Experiment '_patternMode_' cfg.modeMetric '_1_' cfg.hemisphere '_' method '.mat']);

end

function gradientMesh = load_gradient_mesh(cfg, geometry)

if ~exist(cfg.gradientFile, 'file')
    error('Gradient file not found: %s', cfg.gradientFile);
end

S = load(cfg.gradientFile, 'grad_lh_mesh', 'grad_rh_mesh');
gradientMesh = [conv2(S.grad_lh_mesh, ones(2) ./ 4, 'valid') ...
    conv2(S.grad_rh_mesh, ones(2) ./ 4, 'valid')];

if ~isequal(size(gradientMesh), size(geometry.V_mask))
    error('Gradient mesh size does not match the velocity-field mask.');
end

end

function modeData = load_required_mode(cfg, method)

filePath = mode_file(cfg, method);
if ~exist(filePath, 'file')
    error('Principal mode file not found for %s: %s', method, filePath);
end
modeData = load(filePath, 'u1', 's1', 'v1', 'vi_length', 'subjectList');

end

function [DivMap, CurlMap] = mode_to_maps(u1, geometry)

nMode = numel(u1);
if rem(nMode, 2) ~= 0
    error('Mode vector length must be even because it stores [Div; Curl].');
end
halfPoint = nMode / 2;

DivMap = interpolate_to_grid(u1(1:halfPoint), geometry.x, geometry.y, ...
    geometry.xTarget, geometry.yTarget, geometry.V_mask);
CurlMap = interpolate_to_grid(u1(halfPoint+1:end), geometry.x, geometry.y, ...
    geometry.xTarget, geometry.yTarget, geometry.V_mask);

end

function [DivMap, CurlMap] = zscore_mode_maps(DivMap, CurlMap, V_mask)

DivValues = DivMap(V_mask == 1);
CurlValues = CurlMap(V_mask == 1);
DivMap = (DivMap - mean(DivValues, 'omitnan')) ./ std(DivValues, 0, 'omitnan');
CurlMap = (CurlMap - mean(CurlValues, 'omitnan')) ./ std(CurlValues, 0, 'omitnan');

end

function signedMap = threshold_signed_map(mapIn, V_mask, nStd)

signedMap = mapIn;
signedMap(signedMap < nStd & signedMap > -nStd) = 0;
signedMap(signedMap > 0) = 1;
signedMap(signedMap < 0) = -1;
signedMap(V_mask == 0) = nan;

end

function [binnedValues, binnedGradient] = bin_by_gradient(values, gradientVec, nBins)

if isvector(values)
    values = values(:);
end

[gradientVec, idx] = sort(gradientVec(:), 'descend');
values = values(idx, :);
nValues = numel(gradientVec);
edges = round(linspace(1, nValues + 1, nBins + 1));

binnedValues = nan(nBins, size(values, 2));
binnedGradient = nan(nBins, 1);
for i_bin = 1:nBins
    idxBin = edges(i_bin):edges(i_bin + 1) - 1;
    idxBin(idxBin < 1 | idxBin > nValues) = [];
    binnedValues(i_bin, :) = mean(values(idxBin, :), 1, 'omitnan');
    binnedGradient(i_bin) = mean(gradientVec(idxBin), 'omitnan');
end

end

function timeWindow = valid_time_window(requestedWindow, nTime)

timeWindow = requestedWindow(requestedWindow >= 1 & requestedWindow <= nTime);
if isempty(timeWindow)
    error('Requested time window is outside available time points 1:%d.', nTime);
end

end
