function Outputs = rest_global_states(cfg)
%EXFIG3_DIVERSE_MOTIF_GLOBAL_STATES Render Extended Data Fig. 6 global states.
%
% This entry point clusters framewise COF divergence/curl maps into global
% states only when the requested clustering MAT file is absent. Otherwise it
% loads that file and renders the figure directly. The default result path is
% the historical K=20 file already placed in the project-level group directory:
%   data/results/HCP-YA/rfMRI_REST1_LR/group/
%     task-rfMRI_REST1_LR_patterns_Kmeans_correlation_K=20.mat
%
% Example
%   run('environment/matlab_setup.m')
%   outputs = rest_global_states;
%
% A missing result causes a full cohort computation. It can require substantial
% RAM because correlation k-means operates on every selected frame; it never
% overwrites a pre-existing clustering result.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end
cfg = make_config(cfg);

[mask, oi_x, oi_y, V_mask] = make_lr_mask(cfg.parcellation, cfg.dataRoot);
if exist(cfg.clusteringFile, 'file') == 2
    fprintf('Extended Data Fig. 6 | using existing clustering result: %s\n', ...
        cfg.clusteringFile);
    S = load_clustering_result(cfg.clusteringFile, cfg.K, nnz(V_mask));
    computed = false;
else
    fprintf('Extended Data Fig. 6 | clustering result is absent; computing K=%d.\n', cfg.K);
    [S, subjectList] = compute_clustering(cfg, V_mask);
    save_clustering_result(cfg.clusteringFile, S.idx, S.C, subjectList, cfg);
    computed = true;
end

Outputs = struct();
Outputs.cfg = cfg;
Outputs.clusteringFile = cfg.clusteringFile;
Outputs.computedClustering = computed;
Outputs.idx = S.idx;
Outputs.C = S.C;

if cfg.makeFigures
    maps=struct('Div',{},'Curl',{},'Ux',{},'Uy',{});
    nVertices=nnz(V_mask);
    for state=1:cfg.K
        DivMap=zscore_map(vector_to_grid(S.C(state,1:nVertices).',V_mask),V_mask);
        CurlMap=zscore_map(vector_to_grid(S.C(state,nVertices+1:end).',V_mask),V_mask);
        [Ux,Uy]=reconstruct_velocity(DivMap,CurlMap,V_mask,1e-3,1.1);
        maps(state)=struct('Div',DivMap,'Curl',CurlMap,'Ux',Ux,'Uy',Uy);
    end
    Outputs.renderFile=write_plot_data('exFig3_diverse_motif_global_states_render_saved', {maps,mask,oi_x,oi_y,V_mask,cfg,fullfile(cfg.figureDir,cfg.figureFileStem)},7);
else
    Outputs.figure = [];
    Outputs.figureFiles = {};
end
end

function cfg = make_config(overrides)

cfg = rest_config(overrides);
cfg.K = get_default(overrides, 'K', 20);
cfg.startFrame = get_default(overrides, 'startFrame', 101);
cfg.frameCount = get_default(overrides, 'frameCount', 300);
cfg.method = get_default(overrides, 'method', 'COF');
cfg.filtSigma = get_default(overrides, 'filtSigma', 3);
cfg.zThreshold = get_default(overrides, 'zThreshold', 2);
cfg.randomSeed = get_default(overrides, 'randomSeed', 0);
cfg.kmeansReplicates = get_default(overrides, 'kmeansReplicates', 1);
cfg.makeFigures = get_default(overrides, 'makeFigures', true);
cfg.saveFigures = get_default(overrides, 'saveFigures', true);
cfg.closeFigures = get_default(overrides, 'closeFigures', false);
cfg.figureFormats = get_default(overrides, 'figureFormats', {'png', 'pdf'});
cfg.figureFileStem = get_default(overrides, 'figureFileStem', ...
    'rest_global_states');
cfg.clusteringFile = get_default(overrides, 'clusteringFile', fullfile( ...
    cfg.groupMetricDir, sprintf('task-%s_patterns_Kmeans_correlation_K=%d.mat', ...
    cfg.Experiment, cfg.K)));
cfg.leftVelocityFile = get_default(overrides, 'leftVelocityFile', ...
    ['VelField_lh_' cfg.method '.mat']);
cfg.rightVelocityFile = get_default(overrides, 'rightVelocityFile', ...
    ['VelField_rh_' cfg.method '.mat']);
end

function value = get_default(overrides, fieldName, defaultValue)
if isfield(overrides, fieldName) && ~isempty(overrides.(fieldName))
    value = overrides.(fieldName);
else
    value = defaultValue;
end
end

function [mask, oi_x, oi_y, V_mask] = make_lr_mask(parcellation, dataRoot)
[mask_L, ~, oi_x_L, oi_y_L, V_mask_L] = load_cortical_atlas(parcellation, 'L', dataRoot);
[mask_R, ~, oi_x_R, oi_y_R, V_mask_R] = load_cortical_atlas(parcellation, 'R', dataRoot);
mask = [mask_L mask_R];
oi_x = [oi_x_L; oi_x_R + size(mask_L, 2) + 1];
oi_y = [oi_y_L; oi_y_R];
V_mask = [V_mask_L V_mask_R];
end

function S = load_clustering_result(fileName, K, nVertices)
requiredVariables = {'idx', 'C'};
availableVariables = who('-file', fileName);
if ~all(ismember(requiredVariables, availableVariables))
    error('Clustering result must contain idx and C: %s', fileName);
end
S = load(fileName, requiredVariables{:});
if size(S.C, 1) ~= K || size(S.C, 2) ~= 2 * nVertices
    error(['Clustering-centre size %s is incompatible with K=%d and the current ' ...
        'cortical mask (%d vertices): %s'], mat2str(size(S.C)), K, nVertices, fileName);
end
if ~isvector(S.idx) || any(S.idx < 1) || any(S.idx > K)
    error('idx must be a vector of state labels in 1:%d: %s', K, fileName);
end
S.idx = S.idx(:);
end

function [S, subjectList] = compute_clustering(cfg, V_mask)
subjectList = list_complete_subjects(cfg);
if isempty(subjectList)
    error('No subjects with both %s and %s were found under: %s', ...
        cfg.leftVelocityFile, cfg.rightVelocityFile, cfg.vvfDataDir);
end

nSubjects = numel(subjectList);
nFrames = nSubjects * cfg.frameCount;
nVertices = nnz(V_mask);
Feature = zeros(nFrames, 2 * nVertices, 'single');
fprintf('Extended Data Fig. 6 | assembling %d frames from %d subjects.\n', ...
    nFrames, nSubjects);

for i_subject = 1:nSubjects
    subjectID = subjectList{i_subject};
    fprintf('Extended Data Fig. 6 | subject %d/%d: %s\n', ...
        i_subject, nSubjects, subjectID);
    leftFile = fullfile(cfg.vvfDataDir, subjectID, cfg.leftVelocityFile);
    rightFile = fullfile(cfg.vvfDataDir, subjectID, cfg.rightVelocityFile);
    L = load(leftFile, 'VelField');
    R = load(rightFile, 'VelField');
    validate_velocity_field(L, leftFile, cfg);
    validate_velocity_field(R, rightFile, cfg);

    frameRange = cfg.startFrame:(cfg.startFrame + cfg.frameCount - 1);
    if frameRange(end) > size(L.VelField.Ux, 3) || frameRange(end) > size(R.VelField.Ux, 3)
        error('Requested frames %d:%d exceed available data for subject %s.', ...
            frameRange(1), frameRange(end), subjectID);
    end
    Ux = cat(2, L.VelField.Ux(:, :, frameRange), R.VelField.Ux(:, :, frameRange));
    Uy = cat(2, L.VelField.Uy(:, :, frameRange), R.VelField.Uy(:, :, frameRange));
    [DivFrames, CurlFrames] = make_state_features(Ux, Uy, V_mask, cfg);
    rows = (i_subject - 1) * cfg.frameCount + (1:cfg.frameCount);
    Feature(rows, :) = [DivFrames, CurlFrames];
end

if ~isempty(cfg.randomSeed)
    rng(cfg.randomSeed, 'twister');
end
[idx, C] = kmeans(Feature, cfg.K, 'Distance', 'correlation', ...
    'Replicates', cfg.kmeansReplicates);
S = struct('idx', idx, 'C', C);
end

function subjectList = list_complete_subjects(cfg)
entries = dir(cfg.vvfDataDir);
isSubject = [entries.isdir] & ~ismember({entries.name}, {'.', '..'});
candidateSubjects = sort({entries(isSubject).name});
hasFields = cellfun(@(subjectID) exist(fullfile(cfg.vvfDataDir, subjectID, ...
    cfg.leftVelocityFile), 'file') == 2 && exist(fullfile(cfg.vvfDataDir, ...
    subjectID, cfg.rightVelocityFile), 'file') == 2, candidateSubjects);
subjectList = candidateSubjects(hasFields);
end

function validate_velocity_field(S, fileName, cfg)
if ~isfield(S, 'VelField') || ~isfield(S.VelField, 'Ux') || ~isfield(S.VelField, 'Uy')
    error('VelField.Ux and VelField.Uy are required: %s', fileName);
end
if size(S.VelField.Ux, 3) < cfg.startFrame + cfg.frameCount - 1 || ...
        size(S.VelField.Uy, 3) < cfg.startFrame + cfg.frameCount - 1
    error('Velocity field has fewer than %d requested frames: %s', ...
        cfg.startFrame + cfg.frameCount - 1, fileName);
end
end

function [DivFrames, CurlFrames] = make_state_features(Ux, Uy, V_mask, cfg)
imageFilter = fspecial('gaussian', round(3 * cfg.filtSigma), cfg.filtSigma);
nVertices = nnz(V_mask);
DivFrames = zeros(cfg.frameCount, nVertices, 'single');
CurlFrames = zeros(cfg.frameCount, nVertices, 'single');

for t = 1:cfg.frameCount
    divFrame = divergence(Ux(:, :, t), Uy(:, :, t));
    divFrame(V_mask == 0) = nan;
    divFrame = nanconv(divFrame, imageFilter);

    speed = hypot(Ux(:, :, t), Uy(:, :, t));
    speed(speed == 0) = nan;
    [curlFrame, ~] = curl(Ux(:, :, t) ./ speed, Uy(:, :, t) ./ speed);
    curlFrame(V_mask == 0) = nan;
    curlFrame = nanconv(curlFrame, imageFilter);

    DivFrames(t, :) = single(divFrame(V_mask == 1));
    CurlFrames(t, :) = single(curlFrame(V_mask == 1));
end

DivFrames = binarize_frames(DivFrames, cfg.zThreshold);
CurlFrames = binarize_frames(CurlFrames, cfg.zThreshold);
end

function frames = binarize_frames(frames, zThreshold)
validValues = frames(~isnan(frames));
if isempty(validValues)
    error('No valid cortical values are available for motif-state clustering.');
end
mu = mean(validValues, 'omitnan');
sigma = std(validValues, 0, 'omitnan');
if sigma == 0 || ~isfinite(sigma)
    error('Cannot z-score a constant or invalid frame stack.');
end
frames = (frames - mu) ./ sigma;
frames(abs(frames) < zThreshold) = 0;
frames(frames >= zThreshold) = 1;
frames(frames <= -zThreshold) = -1;
frames(isnan(frames)) = 0;
end

function save_clustering_result(fileName, idx, C, subjectList, cfg)
parentDir = fileparts(fileName);
if ~exist(parentDir, 'dir')
    mkdir(parentDir);
end
if exist(fileName, 'file') == 2
    error('Refusing to overwrite existing clustering result: %s', fileName);
end
save(fileName, 'idx', 'C', 'subjectList', 'cfg', '-v7.3');
fprintf('Extended Data Fig. 6 | saved clustering result: %s\n', fileName);
end

function Map = zscore_map(Map, V_mask)
Map(V_mask == 0) = nan;
values = Map(V_mask == 1);
Map = (Map - mean(values, 'omitnan')) ./ std(values, 0, 'omitnan');
end
