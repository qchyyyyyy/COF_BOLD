function rest_surrogate_metrics(i_batch, batch_size, cfg)
%STAT_SURROGATE_MOTIFS_HCP_YA Validate HCP-YA motifs against surrogate VVFs.
%
% This analysis starts from existing velocity fields:
%   data/HCP-YA/rfMRI_REST1_LR/<sub_id>/VelField_lh_COF.mat
%   data/HCP-YA/rfMRI_REST1_LR/<sub_id>/VelField_rh_COF.mat
%   data/HCP-YA/rfMRI_REST1_LR/<sub_id>/surrogate/VelField_lh_COF.mat
%   data/HCP-YA/rfMRI_REST1_LR/<sub_id>/surrogate/VelField_rh_COF.mat
%
% It does not regenerate surrogate data or velocity fields.
%
% Usage:
%   rest_surrogate_metrics()
%   rest_surrogate_metrics(cfg)
%   rest_surrogate_metrics(i_batch, batch_size)
%   rest_surrogate_metrics(i_batch, batch_size, cfg)

if nargin == 1 && isstruct(i_batch)
    cfg = i_batch;
    i_batch = [];
    batch_size = [];
elseif nargin < 3 || isempty(cfg)
    cfg = struct();
end
if nargin < 1
    i_batch = [];
end
if nargin < 2
    batch_size = [];
end
isBatchRun = ~isempty(i_batch);
userSetAggregateStats = isfield(cfg, 'aggregateStats') && ~isempty(cfg.aggregateStats);

repoRoot = find_repo_root(fileparts(mfilename('fullpath')));
cfg = set_default_cfg(cfg, repoRoot);

if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end
if ~exist(cfg.subjectMetricDir, 'dir')
    mkdir(cfg.subjectMetricDir);
end

templateFile = fullfile(cfg.dataRoot, 'templates', 'coordinate_cii_v2.mat');
if ~exist(templateFile, 'file')
    error('Template file not found: %s', templateFile);
end
load(templateFile, 'cii_mask_left', 'cii_mask_right');
[V_mask_L, V_mask_R, V_mask] = make_v_masks(cii_mask_left, cii_mask_right);

dataDir = fullfile(cfg.dataRoot, cfg.Site, cfg.Experiment);
if ~exist(dataDir, 'dir')
    error('Data directory not found: %s', dataDir);
end

subjectList = list_subjects(dataDir);
if isempty(subjectList)
    error('No subject directories found under: %s', dataDir);
end

idx_sub = select_subject_indices(i_batch, batch_size, numel(subjectList));
fprintf('Processing %d/%d subjects in this batch.\n', numel(idx_sub), numel(subjectList));

for i_sub = idx_sub
    subjectID = subjectList{i_sub};
    subjectDir = fullfile(dataDir, subjectID);
    surrogateDir = fullfile(subjectDir, 'surrogate');
    subjectMetricFile = subject_metric_mat_file(cfg, subjectID);
    subjectMetricCsv = subject_metric_csv_file(cfg, subjectID);

    fprintf('Subject %d/%d: %s\n', i_sub, numel(subjectList), subjectID);

    if exist(subjectMetricFile, 'file') && ~cfg.forceSubjectMetrics
        fprintf('Subject metrics already exist: %s\n', subjectMetricFile);
        continue
    end

    realFiles = make_vvf_files(subjectDir, cfg.Method);
    surrogateFiles = make_vvf_files(surrogateDir, cfg.Method);

    missingRows = init_missing_rows();
    [isComplete, missingRows] = check_subject_files(subjectID, realFiles, surrogateFiles, missingRows);
    if ~isComplete
        missingTable = rows_to_missing_table(missingRows);
        writetable(missingTable, fullfile(cfg.subjectMetricDir, [subjectID '_surrogate_motif_missing.csv']));
        continue
    end

    realMetrics = compute_one_dataset(realFiles, V_mask_L, V_mask_R, V_mask, cfg);
    surrogateMetrics = compute_one_dataset(surrogateFiles, V_mask_L, V_mask_R, V_mask, cfg);

    metricRows = init_metric_rows();
    metricRows = append_metric_rows(metricRows, subjectID, 'real', realMetrics);
    metricRows = append_metric_rows(metricRows, subjectID, 'surrogate', surrogateMetrics);

    subjectMetricsTable = rows_to_metric_table(metricRows);
    save(subjectMetricFile, 'subjectMetricsTable', 'subjectID', 'cfg', '-v7.3');
    writetable(subjectMetricsTable, subjectMetricCsv);
end

if isBatchRun && ~userSetAggregateStats
    fprintf('Batch subject metrics saved. Skipping shared aggregation for parallel safety.\n');
    fprintf('Run rest_surrogate_metrics() once after all batches finish to aggregate results.\n');
    return
end

[metricsTable, missingTable] = aggregate_subject_results(subjectList, dataDir, cfg);

writetable(metricsTable, fullfile(cfg.outputDir, 'surrogate_motif_subject_metrics.csv'));
writetable(missingTable, fullfile(cfg.outputDir, 'surrogate_motif_missing_subjects.csv'));

if cfg.aggregateStats && height(metricsTable) > 0
    statsTable = compute_paired_stats(metricsTable, cfg);
    writetable(statsTable, fullfile(cfg.outputDir, 'surrogate_motif_paired_stats.csv'));

    save(fullfile(cfg.outputDir, 'surrogate_motif_validation_results.mat'), ...
        'metricsTable', 'missingTable', 'statsTable', 'cfg', '-v7.3');

    if cfg.makeFigure
        %%
        plot_surrogate_summary(metricsTable, statsTable, cfg);
    end
else
    statsTable = table();
    save(fullfile(cfg.outputDir, 'surrogate_motif_validation_results.mat'), ...
        'metricsTable', 'missingTable', 'statsTable', 'cfg', '-v7.3');
end

fprintf('Saved surrogate motif validation outputs to: %s\n', cfg.outputDir);

end

function cfg = set_default_cfg(cfg, repoRoot)

cfg = set_default(cfg, 'dataRoot', fullfile(repoRoot, 'data'));
cfg = set_default(cfg, 'Site', 'HCP-YA');
cfg = set_default(cfg, 'Experiment', 'rfMRI_REST1_LR');
cfg = set_default(cfg, 'Method', 'COF');
cfg = set_default(cfg, 'parcellation', 'HCPex_360');
cfg = set_default(cfg, 'filtSigma', 3);
cfg = set_default(cfg, 'do_zscore', 1);
cfg = set_default(cfg, 'do_binarize', 1);
cfg = set_default(cfg, 'n_std', 2);
cfg = set_default(cfg, 'size_threshold', 10);
cfg = set_default(cfg, 'nPerm', 10000);
cfg = set_default(cfg, 'makeFigure', true);
cfg = set_default(cfg, 'aggregateStats', true);
cfg = set_default(cfg, 'forceSubjectMetrics', false);
cfg = set_default(cfg, 'randomSeed', 20260709);
cfg = set_default(cfg, 'outputDir', fullfile(repoRoot, 'results', ...
    'surrogate_motif_validation', 'HCP-YA', 'rfMRI_REST1_LR'));
cfg = set_default(cfg, 'subjectMetricDir', fullfile(cfg.outputDir, 'subject_metrics'));

end

function cfg = set_default(cfg, fieldName, defaultValue)

if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = defaultValue;
end

end

function [V_mask_L, V_mask_R, V_mask] = make_v_masks(cii_mask_left, cii_mask_right)

V_mask_L = conv2(cii_mask_left, ones(2, 2) ./ 4, 'valid');
V_mask_R = conv2(cii_mask_right, ones(2, 2) ./ 4, 'valid');
V_mask_L(V_mask_L ~= 1) = 0;
V_mask_R(V_mask_R ~= 1) = 0;
V_mask = [V_mask_L V_mask_R];

end

function subjectList = list_subjects(dataDir)

items = dir(dataDir);
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

function idx_sub = select_subject_indices(i_batch, batch_size, nSubjects)

if isempty(i_batch)
    idx_sub = 1:nSubjects;
    return
end

if isempty(batch_size)
    batch_size = 1;
end

if i_batch < 1 || batch_size < 1
    error('i_batch and batch_size must be positive integers.');
end

idx_sub = i_batch:batch_size:nSubjects;

end

function filePath = subject_metric_mat_file(cfg, subjectID)

filePath = fullfile(cfg.subjectMetricDir, [subjectID '_surrogate_motif_metrics.mat']);

end

function filePath = subject_metric_csv_file(cfg, subjectID)

filePath = fullfile(cfg.subjectMetricDir, [subjectID '_surrogate_motif_metrics.csv']);

end

function [metricsTable, missingTable] = aggregate_subject_results(subjectList, dataDir, cfg)

metricTables = {};
missingRows = init_missing_rows();

for i_sub = 1:numel(subjectList)
    subjectID = subjectList{i_sub};
    subjectDir = fullfile(dataDir, subjectID);
    surrogateDir = fullfile(subjectDir, 'surrogate');

    realFiles = make_vvf_files(subjectDir, cfg.Method);
    surrogateFiles = make_vvf_files(surrogateDir, cfg.Method);
    [isComplete, missingRows] = check_subject_files(subjectID, realFiles, surrogateFiles, missingRows);
    if ~isComplete
        continue
    end

    subjectMetricFile = subject_metric_mat_file(cfg, subjectID);
    if ~exist(subjectMetricFile, 'file')
        missingRows(end+1).subject_id = subjectID; %#ok<AGROW>
        missingRows(end).missing_item = 'subject_metric_file';
        missingRows(end).path = subjectMetricFile;
        continue
    end

    S = load(subjectMetricFile, 'subjectMetricsTable');
    if ~isfield(S, 'subjectMetricsTable')
        missingRows(end+1).subject_id = subjectID; %#ok<AGROW>
        missingRows(end).missing_item = 'subject_metric_table';
        missingRows(end).path = subjectMetricFile;
        continue
    end

    metricTables{end+1, 1} = S.subjectMetricsTable; %#ok<AGROW>
end

if isempty(metricTables)
    metricsTable = rows_to_metric_table(init_metric_rows());
else
    metricsTable = vertcat(metricTables{:});
end

missingTable = rows_to_missing_table(missingRows);

end

function files = make_vvf_files(baseDir, method)

files.lh = fullfile(baseDir, ['VelField_lh_' method '.mat']);
files.rh = fullfile(baseDir, ['VelField_rh_' method '.mat']);

end

function [isComplete, missingRows] = check_subject_files(subjectID, realFiles, surrogateFiles, missingRows)

isComplete = true;
checks = {
    'real_lh', realFiles.lh;
    'real_rh', realFiles.rh;
    'surrogate_lh', surrogateFiles.lh;
    'surrogate_rh', surrogateFiles.rh
    };

for i_check = 1:size(checks, 1)
    if ~exist(checks{i_check, 2}, 'file')
        missingRows(end+1).subject_id = subjectID; %#ok<AGROW>
        missingRows(end).missing_item = checks{i_check, 1};
        missingRows(end).path = checks{i_check, 2};
        isComplete = false;
    end
end

end

function metrics = compute_one_dataset(files, V_mask_L, V_mask_R, V_mask, cfg)

L = load(files.lh, 'VelField');
R = load(files.rh, 'VelField');

Ux = cat(2, L.VelField.Ux, R.VelField.Ux);
Uy = cat(2, L.VelField.Uy, R.VelField.Uy);
clearvars L R

if ~isequal(size(Ux, 1), size(V_mask, 1)) || ~isequal(size(Ux, 2), size(V_mask, 2))
    error('Combined VVF size [%d %d] does not match V_mask size [%d %d].', ...
        size(Ux, 1), size(Ux, 2), size(V_mask, 1), size(V_mask, 2));
end

[DivMap, CurlMap] = velocity_derivatives(Ux, Uy, V_mask, cfg.filtSigma, ...
    cfg.do_zscore, cfg.do_binarize, cfg.n_std);
clearvars Ux Uy

metrics = init_metrics();
metrics = fill_pattern_metrics(metrics, DivMap, CurlMap, V_mask, cfg);

% Keep the per-hemisphere masks referenced to avoid accidental interface drift.
if isempty(V_mask_L) || isempty(V_mask_R)
    error('Invalid hemisphere masks');
end

end

function metrics = init_metrics()

motifs = motif_names();
metricNames = metric_names();
for i_motif = 1:numel(motifs)
    for i_metric = 1:numel(metricNames)
        metrics.(motifs{i_motif}).(metricNames{i_metric}) = nan;
    end
end

end

function metrics = fill_pattern_metrics(metrics, DivMap, CurlMap, V_mask, cfg)

T = size(DivMap, 3);
NumSource = zeros(T, 1);
NumSink = zeros(T, 1);
NumCCSpiral = zeros(T, 1);
NumCSpiral = zeros(T, 1);

SizeSource = zeros(T, 1);
SizeSink = zeros(T, 1);
SizeCCSpiral = zeros(T, 1);
SizeCSpiral = zeros(T, 1);

for t = 1:T
    [NumSource(t), NumSink(t)] = count_motifs(DivMap(:, :, t), cfg.size_threshold);
    [NumCCSpiral(t), NumCSpiral(t)] = count_motifs(CurlMap(:, :, t), cfg.size_threshold);

    [SizeSource(t), SizeSink(t)] = measure_motif_sizes(DivMap(:, :, t));
    [SizeCCSpiral(t), SizeCSpiral(t)] = measure_motif_sizes(CurlMap(:, :, t));
end

metrics.source.PQ = mean(NumSource, 'omitnan');
metrics.sink.PQ = mean(NumSink, 'omitnan');
metrics.CC_spiral.PQ = mean(NumCCSpiral, 'omitnan');
metrics.C_spiral.PQ = mean(NumCSpiral, 'omitnan');

metrics.source.PS = mean(SizeSource, 'omitnan');
metrics.sink.PS = mean(SizeSink, 'omitnan');
metrics.CC_spiral.PS = mean(SizeCCSpiral, 'omitnan');
metrics.C_spiral.PS = mean(SizeCSpiral, 'omitnan');

metrics.source.dwell_time = mean_dwell_time(DivMap, 1, cfg.size_threshold);
metrics.sink.dwell_time = mean_dwell_time(DivMap, -1, cfg.size_threshold);
metrics.CC_spiral.dwell_time = mean_dwell_time(CurlMap, 1, cfg.size_threshold);
metrics.C_spiral.dwell_time = mean_dwell_time(CurlMap, -1, cfg.size_threshold);

metrics.source.spatial_sd = motif_spatial_sd(DivMap == 1, V_mask);
metrics.sink.spatial_sd = motif_spatial_sd(DivMap == -1, V_mask);
metrics.CC_spiral.spatial_sd = motif_spatial_sd(CurlMap == 1, V_mask);
metrics.C_spiral.spatial_sd = motif_spatial_sd(CurlMap == -1, V_mask);

end

function value = mean_dwell_time(Map, motifValue, sizeThreshold)

[~, DwellTime] = track_motifs(Map, motifValue, sizeThreshold);
if isempty(DwellTime)
    value = 0;
else
    value = mean(DwellTime, 'omitnan');
end

end

function value = motif_spatial_sd(MotifMap, V_mask)

T = size(MotifMap, 3);
MotifMap(repmat(V_mask, 1, 1, T) == 0) = false;
densityMap = mean(double(MotifMap), 3, 'omitnan');
value = std(densityMap(V_mask == 1), [], 'omitnan');

end

function metricRows = append_metric_rows(metricRows, subjectID, dataType, metrics)

motifs = motif_names();
metricNames = metric_names();

for i_motif = 1:numel(motifs)
    motif = motifs{i_motif};
    for i_metric = 1:numel(metricNames)
        metric = metricNames{i_metric};
        metricRows(end+1).subject_id = subjectID; %#ok<AGROW>
        metricRows(end).data_type = dataType;
        metricRows(end).motif = motif;
        metricRows(end).metric = metric;
        metricRows(end).value = metrics.(motif).(metric);
    end
end

end

function metricRows = init_metric_rows()

metricRows = struct('subject_id', {}, 'data_type', {}, 'motif', {}, 'metric', {}, 'value', {});

end

function missingRows = init_missing_rows()

missingRows = struct('subject_id', {}, 'missing_item', {}, 'path', {});

end

function T = rows_to_metric_table(rows)

if isempty(rows)
    T = table(cell(0, 1), cell(0, 1), cell(0, 1), cell(0, 1), zeros(0, 1), ...
        'VariableNames', {'subject_id', 'data_type', 'motif', 'metric', 'value'});
else
    T = struct2table(rows);
end

end

function T = rows_to_missing_table(rows)

if isempty(rows)
    T = table(cell(0, 1), cell(0, 1), cell(0, 1), ...
        'VariableNames', {'subject_id', 'missing_item', 'path'});
else
    T = struct2table(rows);
end

end

function statsTable = compute_paired_stats(metricsTable, cfg)

motifs = motif_names();
metricNames = metric_names();
statsRows = struct('motif', {}, 'metric', {}, 'n', {}, 'real_mean', {}, ...
    'surrogate_mean', {}, 'mean_diff', {}, 'cohen_dz', {}, 'p_perm', {}, 'q_fdr', {});

rng(cfg.randomSeed);

for i_motif = 1:numel(motifs)
    for i_metric = 1:numel(metricNames)
        [realVals, surrogateVals] = paired_values(metricsTable, motifs{i_motif}, metricNames{i_metric});
        [p_perm, meanDiff, cohenDz, nPair] = paired_signflip_test(realVals, surrogateVals, cfg.nPerm);

        statsRows(end+1).motif = motifs{i_motif}; %#ok<AGROW>
        statsRows(end).metric = metricNames{i_metric};
        statsRows(end).n = nPair;
        statsRows(end).real_mean = mean(realVals, 'omitnan');
        statsRows(end).surrogate_mean = mean(surrogateVals, 'omitnan');
        statsRows(end).mean_diff = meanDiff;
        statsRows(end).cohen_dz = cohenDz;
        statsRows(end).p_perm = p_perm;
        statsRows(end).q_fdr = nan;
    end
end

statsTable = struct2table(statsRows);
statsTable.q_fdr = bh_fdr(statsTable.p_perm);

end

function [realVals, surrogateVals] = paired_values(metricsTable, motif, metric)

subjectIDs = unique(metricsTable.subject_id, 'stable');
realVals = nan(numel(subjectIDs), 1);
surrogateVals = nan(numel(subjectIDs), 1);

for i_sub = 1:numel(subjectIDs)
    subID = subjectIDs{i_sub};
    idxReal = strcmp(metricsTable.subject_id, subID) & ...
        strcmp(metricsTable.data_type, 'real') & ...
        strcmp(metricsTable.motif, motif) & ...
        strcmp(metricsTable.metric, metric);
    idxSurrogate = strcmp(metricsTable.subject_id, subID) & ...
        strcmp(metricsTable.data_type, 'surrogate') & ...
        strcmp(metricsTable.motif, motif) & ...
        strcmp(metricsTable.metric, metric);

    if any(idxReal)
        realVals(i_sub) = metricsTable.value(find(idxReal, 1, 'first'));
    end
    if any(idxSurrogate)
        surrogateVals(i_sub) = metricsTable.value(find(idxSurrogate, 1, 'first'));
    end
end

valid = ~isnan(realVals) & ~isnan(surrogateVals);
realVals = realVals(valid);
surrogateVals = surrogateVals(valid);

end

function [p_perm, meanDiff, cohenDz, nPair] = paired_signflip_test(realVals, surrogateVals, nPerm)

diffVals = realVals - surrogateVals;
diffVals = diffVals(~isnan(diffVals));
nPair = numel(diffVals);

if nPair == 0
    p_perm = nan;
    meanDiff = nan;
    cohenDz = nan;
    return
end

meanDiff = mean(diffVals, 'omitnan');
diffStd = std(diffVals, [], 'omitnan');
if nPair < 2 || diffStd == 0
    cohenDz = nan;
else
    cohenDz = meanDiff ./ diffStd;
end

nullMeans = zeros(nPerm, 1);
for i_perm = 1:nPerm
    signs = double(rand(nPair, 1) > 0.5);
    signs(signs == 0) = -1;
    nullMeans(i_perm) = mean(diffVals .* signs, 'omitnan');
end

p_perm = (sum(abs(nullMeans) >= abs(meanDiff)) + 1) ./ (nPerm + 1);

end

function qVals = bh_fdr(pVals)

qVals = nan(size(pVals));
valid = ~isnan(pVals);
p = pVals(valid);

if isempty(p)
    return
end

[pSorted, order] = sort(p);
m = numel(pSorted);
qSorted = pSorted .* m ./ (1:m)';
qSorted = flipud(cummin(flipud(qSorted)));
qSorted(qSorted > 1) = 1;

q = nan(size(p));
q(order) = qSorted;
qVals(valid) = q;

end

function plot_surrogate_summary(metricsTable, statsTable, cfg)

metricNames = metric_names();
motifs = motif_names();
motifLabels = {'Diverging', 'Converging', 'CC-spiral', 'C-spiral'};
realColor = [0.22 0.54 0.70];
surrogateColor = [0.70 0.70 0.70];
hasBoxchart = exist('boxchart', 'file') == 2 || exist('boxchart', 'builtin') == 5;
hasTiledlayout = exist('tiledlayout', 'file') == 2 || exist('tiledlayout', 'builtin') == 5;

fig = figure('Color', 'w', 'Position', [100 100 1200 760]);
if hasTiledlayout
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
end

for i_metric = 1:numel(metricNames)
    metric = metricNames{i_metric};
    if hasTiledlayout
        nexttile;
    else
        subplot(2, 2, i_metric);
    end

    qVals = nan(numel(motifs), 1);
    allVals = [];
    hReal = [];
    hSurrogate = [];

    for i_motif = 1:numel(motifs)
        [realVals, surrogateVals] = paired_values(metricsTable, motifs{i_motif}, metric);
        allVals = [allVals; realVals(:); surrogateVals(:)]; %#ok<AGROW>

        xReal = (i_motif - 0.18) .* ones(size(realVals));
        xSurrogate = (i_motif + 0.18) .* ones(size(surrogateVals));
        if hasBoxchart
            hRealTmp = boxchart(xReal, realVals, ...
                'BoxFaceColor', realColor, ...
                'MarkerColor', realColor, ...
                'BoxWidth', 0.28);
            hold on
            hSurrogateTmp = boxchart(xSurrogate, surrogateVals, ...
                'BoxFaceColor', surrogateColor, ...
                'MarkerColor', [0.45 0.45 0.45], ...
                'BoxWidth', 0.28);
        else
            boxplot(realVals, 'positions', i_motif - 0.18, 'widths', 0.28, ...
                'colors', realColor, 'symbol', '');
            hold on
            boxplot(surrogateVals, 'positions', i_motif + 0.18, 'widths', 0.28, ...
                'colors', surrogateColor, 'symbol', '');
            hRealTmp = plot(nan, nan, 's', 'MarkerFaceColor', realColor, ...
                'MarkerEdgeColor', realColor);
            hSurrogateTmp = plot(nan, nan, 's', 'MarkerFaceColor', surrogateColor, ...
                'MarkerEdgeColor', surrogateColor);
        end
        if isempty(hReal)
            hReal = hRealTmp;
            hSurrogate = hSurrogateTmp;
        end

        idxStat = strcmp(statsTable.motif, motifs{i_motif}) & strcmp(statsTable.metric, metric);
        if any(idxStat)
            qVals(i_motif) = statsTable.q_fdr(find(idxStat, 1, 'first'));
        end
    end

    [yMin, yMax, yPad] = data_axis_limits(allVals);
    for i_motif = 1:numel(motifs)
        text(i_motif, yMax + 0.45 .* yPad, ...
            q_label(qVals(i_motif)), 'HorizontalAlignment', 'center', 'FontSize', 8);
    end

    set(gca, 'XTick', 1:numel(motifs), 'XTickLabel', motifLabels, ...
        'XTickLabelRotation', 25, 'LineWidth', 1, 'FontSize', 9);
    xlim([0.45 numel(motifs) + 0.55]);
    ylabel(metric, 'Interpreter', 'none');
    ylim([yMin - yPad, yMax + yPad]);
    if i_metric == 1
        legend([hReal hSurrogate], {'Real', 'Surrogate'}, 'Location', 'best', 'Box', 'off');
    end
    box off
end

figurePathPng = fullfile(cfg.outputDir, 'surrogate_motif_validation_summary.png');
figurePathPdf = fullfile(cfg.outputDir, 'surrogate_motif_validation_summary.pdf');

try
    exportgraphics(fig, figurePathPng, 'Resolution', 300);
    exportgraphics(fig, figurePathPdf, 'ContentType', 'vector');
catch
    saveas(fig, figurePathPng);
    saveas(fig, figurePathPdf);
end

% close(fig);

end

function [yMin, yMax, yPad] = data_axis_limits(values)

values = values(~isnan(values));
if isempty(values)
    yMin = 0;
    yMax = 1;
    yPad = 0.1;
else
    yMin = min(values);
    yMax = max(values);
    yRange = yMax - yMin;
    if yRange == 0
        yRange = max(abs(yMax), 1);
    end
    yPad = 0.12 .* yRange;
end

end

function label = q_label(q)

if isnan(q)
    label = 'q=n/a';
elseif q < 0.001
    label = 'q<0.001';
else
    label = sprintf('q=%.3f', q);
end

end

function motifs = motif_names()

motifs = {'source', 'sink', 'CC_spiral', 'C_spiral'};

end

function metricNames = metric_names()

metricNames = {'PQ', 'PS', 'dwell_time', 'spatial_sd'};

end

function repoRoot = find_repo_root(startDir)

repoRoot = startDir;
while true
    if isfile(fullfile(repoRoot, 'setup.m'))
        return
    end
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error('Could not locate repository root from: %s', startDir);
    end
    repoRoot = parentDir;
end

end
