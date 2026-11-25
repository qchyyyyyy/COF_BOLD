function Outputs = rest_motif_densities(cfg)
%FIG1E_STATISTIC_PATTERNSP Summarize Fig. 1e spatial motif-density maps.
%
% The input priority is:
%   1) an explicitly configured compact MAT file containing meanMaps;
%   2) group-level numPattern_voxel files produced by rest_motif_metrics with
%      save_method='group';
%   3) participant-level numPattern_voxel files as a fallback.
%
% The default data root is resolved by rest_config. During the current
% project phase it is <program-root>/data; the later reproduction-local data
% migration is controlled at that single shared configuration switch point.

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

cfg = rest_config(cfg);

[~, ~, oi_x, oi_y, V_mask] = make_lr_mask(cfg.parcellation, cfg.dataRoot);
nMethods = numel(cfg.Method);
nPatterns = numel(cfg.patternSourceNames);
[meanMaps, patternMaps, subjectList, usedAggregatedMaps] = load_or_compute_mean_maps(cfg, V_mask, nMethods, nPatterns);

Outputs = struct();
Outputs.cfg = cfg;
Outputs.subjectList = subjectList;
Outputs.V_mask = V_mask;
Outputs.patternMaps = patternMaps;
Outputs.meanMaps = meanMaps;

if cfg.saveData
    if ~exist(cfg.figureDataDir, 'dir')
        mkdir(cfg.figureDataDir);
    end
    outputFile = fullfile(cfg.figureDataDir, 'fig1e_pattern_spatial_maps.mat');
    if usedAggregatedMaps
        save(outputFile, 'meanMaps', 'subjectList', 'cfg', '-v7.3');
    else
        save(outputFile, 'meanMaps', 'patternMaps', 'subjectList', 'cfg', '-v7.3');
    end
end

if cfg.makeFigures
    if ~exist(cfg.figureDir, 'dir')
        mkdir(cfg.figureDir);
    end
    Outputs.renderFile=write_plot_data('fig1e_statistic_patternSP_render_saved', {meanMaps,V_mask,oi_x,oi_y,cfg,fullfile(cfg.figureDir,'fig1_density')},6);
end

end

function [meanMaps, patternMaps, subjectList, usedAggregatedMaps] = load_or_compute_mean_maps(cfg, V_mask, nMethods, nPatterns)

usedAggregatedMaps = isfield(cfg, 'precomputedMeanMapsFile') && ...
    ~isempty(cfg.precomputedMeanMapsFile);
patternMaps = [];

if usedAggregatedMaps
    inputFile = cfg.precomputedMeanMapsFile;
    if ~exist(inputFile, 'file')
        error('Precomputed Fig. 1e mean-map file not found: %s', inputFile);
    end
    S = load(inputFile, 'meanMaps');
    if ~isfield(S, 'meanMaps')
        error('Precomputed Fig. 1e file does not contain meanMaps: %s', inputFile);
    end
    meanMaps = S.meanMaps;
    subjectList = {};
    inputVariables = who('-file', inputFile);
    if ismember('subjectList', inputVariables)
        subjectData = load(inputFile, 'subjectList');
        subjectList = subjectData.subjectList;
    end
    expectedSize = [size(V_mask) nMethods nPatterns];
    if ~isequal(size(meanMaps), expectedSize)
        error('Precomputed meanMaps size %s does not match expected size %s.', ...
            mat2str(size(meanMaps)), mat2str(expectedSize));
    end
    fprintf('Fig. 1e | using precomputed group mean maps: %s\n', inputFile);
    return
end

% A group-level rest_motif_metrics run stores one file per method, with the
% participant dimension in the third array dimension. Prefer these files
% when available so Fig. 1e does not require participant-level output
% directories to be present locally.
[groupFiles, groupPresent, groupMetricDir] = get_group_metric_files(cfg, nMethods);
if any(groupPresent)
    if ~all(groupPresent)
        missingFiles = groupFiles(~groupPresent);
        error(['Group-level numPattern_voxel output is incomplete. Missing: %s. ' ...
            'Complete the group run or configure cfg.precomputedMeanMapsFile.'], ...
            strjoin(missingFiles, '; '));
    end

    [meanMaps, subjectList] = load_group_mean_maps( ...
        cfg, V_mask, nMethods, nPatterns, groupFiles, groupMetricDir);
    usedAggregatedMaps = true;
    fprintf('Fig. 1e | using group-level metric files under: %s\n', ...
        groupMetricDir);
    return
end

subjectList = list_subjects(cfg.participantMetricDir);
if isempty(subjectList)
    error('No participant metric subject directories found under: %s', cfg.participantMetricDir);
end

nSubjects = numel(subjectList);
mapSize = size(V_mask);
patternMaps = zeros([mapSize nSubjects nMethods nPatterns]);

for i_method = 1:nMethods
    method = cfg.Method{i_method};
    for i_sub = 1:nSubjects
        subjectID = subjectList{i_sub};
        fprintf('Fig. 1e | %s | subject %d/%d: %s\n', method, i_sub, nSubjects, subjectID);
        metricFile = fullfile(cfg.participantMetricDir, subjectID, ...
            [subjectID '_' cfg.Experiment '_numPattern_voxel_' cfg.hemisphere '_' method '.mat']);
        S = load_pattern_metric(metricFile, cfg.patternSourceNames);

        for i_pattern = 1:nPatterns
            patternMaps(:, :, i_sub, i_method, i_pattern) = S.(cfg.patternSourceNames{i_pattern});
        end
    end
end

meanMaps = mean(patternMaps, 3, 'omitnan');
meanMaps = reshape(meanMaps, [mapSize nMethods nPatterns]);

end

function [groupFiles, groupPresent, groupMetricDir] = get_group_metric_files(cfg, nMethods)
% Return the group-level files expected from rest_motif_metrics.
%
% The regenerated pipeline uses cfg.groupMetricDir. The second candidate is
% the legacy/shared result layout used by the downloaded HCP group outputs:
%   <dataRoot>/results/<Site>/group/

candidateDirs = {cfg.groupMetricDir};
legacyGroupDir = fullfile(cfg.dataRoot, 'results', cfg.Site, 'group');
if ~strcmp(legacyGroupDir, cfg.groupMetricDir)
    candidateDirs{end+1} = legacyGroupDir; %#ok<AGROW>
end

partialFiles = {};
partialPresent = false;
for i_dir = 1:numel(candidateDirs)
    candidateFiles = make_group_metric_files(candidateDirs{i_dir}, ...
        cfg, nMethods);
    candidatePresent = cellfun(@(fileName) exist(fileName, 'file') == 2, ...
        candidateFiles);

    if all(candidatePresent)
        groupFiles = candidateFiles;
        groupPresent = candidatePresent;
        groupMetricDir = candidateDirs{i_dir};
        return
    end

    if any(candidatePresent) && ~partialPresent
        partialFiles = candidateFiles;
        partialPresent = candidatePresent;
    end
end

if ~isempty(partialFiles)
    groupFiles = partialFiles;
    groupPresent = partialPresent;
    groupMetricDir = fileparts(groupFiles{1});
else
    groupFiles = make_group_metric_files(cfg.groupMetricDir, cfg, nMethods);
    groupPresent = false(size(groupFiles));
    groupMetricDir = cfg.groupMetricDir;
end
end

function groupFiles = make_group_metric_files(groupMetricDir, cfg, nMethods)

groupFiles = cell(nMethods, 1);
for i_method = 1:nMethods
    groupFiles{i_method} = fullfile(groupMetricDir, ...
        ['task-' cfg.Experiment '_numPattern_voxel_' cfg.hemisphere '_' ...
        cfg.Method{i_method} '.mat']);
end
end

function [meanMaps, subjectList] = load_group_mean_maps( ...
    cfg, V_mask, nMethods, nPatterns, groupFiles, groupMetricDir)
% Load group files and average their participant dimension.

mapSize = size(V_mask);
meanMaps = nan([mapSize nMethods nPatterns]);

for i_method = 1:nMethods
    for i_pattern = 1:nPatterns
        requestedField = cfg.patternSourceNames{i_pattern};
        mapField = resolve_group_map_field(groupFiles{i_method}, ...
            requestedField);
        meanMaps(:, :, i_method, i_pattern) = mean_group_map( ...
            groupFiles{i_method}, mapField, mapSize);
    end
end

subjectList = load_group_subject_list(cfg, groupMetricDir);
end

function mapField = resolve_group_map_field(groupFile, requestedField)
% rest_motif_metrics saves SourceMap/SinkMap/CCVortexMap/CVortexMap.

availableFields = who('-file', groupFile);
if ismember(requestedField, availableFields)
    mapField = requestedField;
    return
end

candidateField = [requestedField 'Map'];
if ismember(candidateField, availableFields)
    mapField = candidateField;
    return
end

error('No map variable %s or %s found in group file: %s', ...
    requestedField, candidateField, groupFile);
end

function mapMean = mean_group_map(groupFile, mapField, mapSize)
% Average a group map without loading all participants into MATLAB memory.

fileInfo = whos('-file', groupFile);
fieldIndex = strcmp({fileInfo.name}, mapField);
if ~any(fieldIndex)
    error('Map variable %s not found in group file: %s', ...
        mapField, groupFile);
end

fieldSize = fileInfo(fieldIndex).size;
if numel(fieldSize) == 2
    loaded = load(groupFile, mapField);
    mapMean = loaded.(mapField);
    if ~isequal(size(mapMean), mapSize)
        error(['Unexpected size for %s in %s. Expected [%s], got [%s].'], ...
            mapField, groupFile, size_string(mapSize), size_string(size(mapMean)));
    end
    return
end

if numel(fieldSize) ~= 3 || ~isequal(fieldSize(1:2), mapSize)
    error(['Unexpected size for %s in %s. Expected [%s] or [%s N_subjects], ' ...
        'got [%s].'], mapField, groupFile, size_string(mapSize), ...
        size_string([mapSize 1]), size_string(fieldSize));
end

groupData = matfile(groupFile);
nSubjects = fieldSize(3);
chunkSize = 16;
sumMap = zeros(mapSize);
countMap = zeros(mapSize);

for firstSubject = 1:chunkSize:nSubjects
    lastSubject = min(firstSubject + chunkSize - 1, nSubjects);
    mapChunk = groupData.(mapField)(1:mapSize(1), 1:mapSize(2), ...
        firstSubject:lastSubject);
    validChunk = isfinite(mapChunk);
    mapChunk(~validChunk) = 0;
    sumMap = sumMap + sum(mapChunk, 3);
    countMap = countMap + sum(validChunk, 3);
end

mapMean = sumMap ./ countMap;
mapMean(countMap == 0) = NaN;
end

function subjectList = load_group_subject_list(cfg, groupMetricDir)
% Load the subject list written by the group aggregation when available.

subjectList = {};
subjectFile = fullfile(groupMetricDir, ...
    [cfg.Experiment '_SubList.mat']);
if exist(subjectFile, 'file') ~= 2
    return
end

availableFields = who('-file', subjectFile);
if ~ismember('subjectList', availableFields)
    return
end

loaded = load(subjectFile, 'subjectList');
subjectList = loaded.subjectList;
end

function text = size_string(value)
text = strtrim(sprintf('%d ', value));
end

function [mask, label, oi_x, oi_y, V_mask] = make_lr_mask(parcellation, dataRoot)

[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = load_cortical_atlas(parcellation, 'L', dataRoot);
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = load_cortical_atlas(parcellation, 'R', dataRoot);

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L; oi_x_R + size(mask_L, 2)];
oi_y = [oi_y_L; oi_y_R];
V_mask = [V_mask_L V_mask_R];

end

function subjectList = list_subjects(participantMetricDir)

if ~exist(participantMetricDir, 'dir')
    subjectList = {};
    return
end

items = dir(participantMetricDir);
items = items([items.isdir]);
subjectList = {};
for i_item = 1:numel(items)
    name = items(i_item).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue
    end
    subjectList{end+1, 1} = name; %#ok<AGROW>
end

subjectList = sort(subjectList);

end

function S = load_pattern_metric(metricFile, patternNames)

if ~exist(metricFile, 'file')
    error('Participant motif-density file not found: %s', metricFile);
end

% rest_motif_metrics stores REST maps with the explicit *Map suffix. Keep the
% figure-level names concise while accepting both legacy and regenerated
% participant files.
S = load(metricFile);
for i_pattern = 1:numel(patternNames)
    patternName = patternNames{i_pattern};
    if isfield(S, patternName)
        continue
    end
    mapName = [patternName 'Map'];
    if isfield(S, mapName)
        S.(patternName) = S.(mapName);
    else
        error('Missing variable %s (or %s) in %s.', patternName, mapName, metricFile);
    end
end

end
