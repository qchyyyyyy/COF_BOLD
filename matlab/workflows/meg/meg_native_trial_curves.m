function summary = meg_native_trial_curves(configFile, subjectId, scanIds)
%MEG_NATIVE_TRIAL_CURVES Build event-aligned trial-time PQ/PS summaries.
%
% This follow-up keeps the P04-EXP-003 scalar PQ/PS interval derivatives
% unchanged, but re-aligns them to every official task image and Fix15sec
% pseudo-trial with a [-0.5, 2.5] s window.  The per-trial table is retained
% on HPC; scan-level time curves are written for participant-aware aggregation.

if nargin < 2 || nargin > 3
    error('MEG:WorkflowArguments', ...
        'Usage: meg_native_trial_curves(configFile, subjectId, [scanIds])');
end
cfg = load_config(configFile);
codeRoot = fileparts(mfilename('fullpath'));
addpath(codeRoot);
addpath(cfg.fieldtrip_root);
ft_defaults;
addpath(fullfile(cfg.megconnectome_root, 'trial_functions'));

subjectId = validate_identifier(subjectId, 'subjectId');
if nargin < 3 || isempty(scanIds)
    scanIds = cfg.scan_ids;
end
scanIds = normalize_scan_ids(scanIds);
if isempty(scanIds)
    error('MEG:WorkflowScans', 'At least one scan must be configured.');
end

resultsRoot = resolve_root(cfg, 'results_root_env', 'results_root');
subjectResultRoot = fullfile(resultsRoot, char(cfg.results_subdir), subjectId);
if ~isfolder(subjectResultRoot), mkdir(subjectResultRoot); end

scanSummary = repmat(struct('scan_id', '', 'status', '', ...
    'raw_file', '', 'metric_file', '', 'trial_file', '', 'scan_time_file', '', ...
    'n_trials', 0, 'n_valid_trials', 0, 'n_time_rows', 0, ...
    'error_identifier', '', 'error_message', ''), 1, numel(scanIds));

for iScan = 1:numel(scanIds)
    scanId = scanIds{iScan};
    scanSummary(iScan).scan_id = scanId;
    scanResultRoot = fullfile(subjectResultRoot, ['scan-' scanId]);
    if ~isfolder(scanResultRoot), mkdir(scanResultRoot); end
    try
        rawFile = resolve_raw_file(cfg, subjectId, scanId);
        metricFile = fullfile(resultsRoot, char(cfg.parent_metric_results_subdir), ...
            subjectId, ['scan-' scanId], 'pq_ps_per_interval.csv');
        sessionDir = resolve_session_directory(cfg, subjectId, scanId);
        metadataFile = fullfile(sessionDir, 'session_metadata.mat');
        required = {rawFile, metricFile, metadataFile, ...
            fullfile(sessionDir, 'session_manifest.json')};
        for iFile = 1:numel(required)
            if ~isfile(required{iFile})
                error('MEG:WorkflowInputMissing', ...
                    'Required input not found: %s', required{iFile});
            end
        end

        loaded = load(metadataFile, 'metadata');
        if ~isfield(loaded, 'metadata') || ...
                ~isfield(loaded.metadata, 'session_input_info')
            error('MEG:WorkflowMetadata', ...
                'session_metadata.mat lacks session_input_info.');
        end
        sessionInfo = loaded.metadata.session_input_info;
        metric = read_metric_table(metricFile, subjectId, scanId);
        decoded = meg_decode_events(rawFile, cfg.megconnectome_root, ...
            cfg.trial_pre_seconds, cfg.trial_post_seconds);
        validate_scan_bounds(decoded, sessionInfo, subjectId, scanId);

        [trialTable, scanTable, qc] = build_trial_time_tables(...
            decoded, metric, sessionInfo, cfg, subjectId, scanId);
        trialFile = fullfile(scanResultRoot, 'trial_time_metrics.csv');
        scanTimeFile = fullfile(scanResultRoot, 'scan_time_means.csv');
        writetable(trialTable, trialFile);
        writetable(scanTable, scanTimeFile);

        alignment = struct;
        alignment.schema_version = 1;
        alignment.status = 'complete';
        alignment.workflow = 'meg_native_trial_curves';
        alignment.parent_metric_experiment = char(cfg.parent_metric_results_subdir);
        alignment.subject_id = subjectId;
        alignment.scan_id = scanId;
        alignment.raw_file = rawFile;
        alignment.metric_file = metricFile;
        alignment.session_metadata_file = metadataFile;
        alignment.recording_fsample_hz = decoded.recording_fsample_hz;
        alignment.recording_samples = decoded.recording_samples;
        alignment.scan_start_sample = decoded.scan_start_sample;
        alignment.scan_end_sample = decoded.scan_end_sample;
        alignment.trial_pre_seconds = cfg.trial_pre_seconds;
        alignment.trial_post_seconds = cfg.trial_post_seconds;
        alignment.frame_interval_seconds = cfg.frame_interval_seconds;
        alignment.n_time_bins = cfg.n_time_bins;
        alignment.n_official_trials = numel(decoded.trials);
        alignment.n_trial_time_rows = height(trialTable);
        alignment.n_scan_time_rows = height(scanTable);
        alignment.qc = qc;
        alignment.decoder_warning_count = decoded.warning_count;
        alignment.created = datestr(now, 30);
        write_json(fullfile(scanResultRoot, 'trial_time_alignment.json'), alignment);
        save(fullfile(scanResultRoot, 'trial_time_alignment.mat'), ...
            'alignment', 'qc');

        scanSummary(iScan).status = 'complete';
        scanSummary(iScan).raw_file = rawFile;
        scanSummary(iScan).metric_file = metricFile;
        scanSummary(iScan).trial_file = trialFile;
        scanSummary(iScan).scan_time_file = scanTimeFile;
        scanSummary(iScan).n_trials = qc.n_trials_total;
        scanSummary(iScan).n_valid_trials = qc.n_valid_trials;
        scanSummary(iScan).n_time_rows = height(scanTable);
    catch ME
        scanSummary(iScan).status = 'failed';
        scanSummary(iScan).error_identifier = ME.identifier;
        scanSummary(iScan).error_message = getReport(ME, 'extended', 'hyperlinks', 'off');
        failure = struct('status', 'failed', 'subject_id', subjectId, ...
            'scan_id', scanId, 'error_identifier', ME.identifier, ...
            'error_message', scanSummary(iScan).error_message, ...
            'created', datestr(now, 30));
        write_json(fullfile(scanResultRoot, 'trial_time_alignment.json'), failure);
    end
end

summary = struct;
summary.schema_version = 1;
summary.workflow = 'meg_native_trial_curves';
summary.subject_id = subjectId;
summary.scan_ids = scanIds;
summary.scans = scanSummary;
summary.status = 'complete';
summary.config_file = configFile;
summary.results_root = resultsRoot;
summary.created = datestr(now, 30);
if any(strcmp({scanSummary.status}, 'failed'))
    summary.status = 'incomplete';
end
write_json(fullfile(subjectResultRoot, 'subject_summary.json'), summary);
save(fullfile(subjectResultRoot, 'subject_summary.mat'), 'summary');
if ~strcmp(summary.status, 'complete')
    error('MEG:WorkflowSubjectFailed', ...
        'meg_native_trial_curves failed for subject %s.', subjectId);
end
end

function cfg = load_config(configFile)
if nargin ~= 1 || ~isfile(configFile)
    error('MEG:WorkflowConfig', 'Configuration file not found: %s', char(configFile));
end
cfg = jsondecode(fileread(configFile));
required = {'workflow','parent_metric_results_subdir', ...
    'raw_root_env','raw_root','processed_root_env','processed_root', ...
    'results_root_env','results_root','fieldtrip_root_env','fieldtrip_root', ...
    'megconnectome_root_env','megconnectome_root','raw_datafile_template', ...
    'results_subdir','scan_ids','frame_interval_seconds','n_time_bins', ...
    'trial_pre_seconds','trial_post_seconds', ...
    'min_trial_coverage_fraction','require_parent_block_containment'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i})
        error('MEG:WorkflowConfigField', 'Missing field %s.', required{i});
    end
end
cfg.raw_root = resolve_root(cfg, 'raw_root_env', 'raw_root');
cfg.processed_root = resolve_root(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_root(cfg, 'results_root_env', 'results_root');
cfg.fieldtrip_root = resolve_root(cfg, 'fieldtrip_root_env', 'fieldtrip_root');
cfg.megconnectome_root = resolve_root(cfg, 'megconnectome_root_env', 'megconnectome_root');
if ~isfield(cfg, 'min_trials_per_timepoint')
    cfg.min_trials_per_timepoint = 1;
end
for f = {'raw_root','processed_root','results_root','fieldtrip_root','megconnectome_root'}
    if ~isfolder(cfg.(f{1}))
        error('MEG:WorkflowRoot', 'Root not found: %s', cfg.(f{1}));
    end
end
expectedBins = round((cfg.trial_pre_seconds + cfg.trial_post_seconds) / ...
    cfg.frame_interval_seconds);
if expectedBins < 1 || cfg.n_time_bins ~= expectedBins
    error('MEG:WorkflowConfigBins', ...
        'n_time_bins must equal round((pre+post)/frame_interval).');
end
if cfg.frame_interval_seconds <= 0 || cfg.trial_pre_seconds < 0 || ...
        cfg.trial_post_seconds <= 0 || cfg.min_trial_coverage_fraction <= 0 || ...
        cfg.min_trial_coverage_fraction > 1
    error('MEG:WorkflowConfigValue', 'Invalid time or coverage parameters.');
end
end

function metric = read_metric_table(filename, subjectId, scanId)
metric = readtable(filename);
required = {'interval_frame_index','interval_time_seconds','PQ','PS'};
for i = 1:numel(required)
    if ~ismember(required{i}, metric.Properties.VariableNames)
        error('MEG:WorkflowMetricColumns', ...
            'Missing metric column %s in %s.', required{i}, filename);
    end
end
if ismember('subject_id', metric.Properties.VariableNames)
    ids = identifiers_to_cellstr(metric.subject_id);
    if any(~strcmp(ids, subjectId))
        error('MEG:WorkflowMetricSubject', 'Metric subject mismatch.');
    end
end
if ismember('scan_id', metric.Properties.VariableNames)
    ids = identifiers_to_cellstr(metric.scan_id);
    if any(~strcmp(ids, scanId))
        error('MEG:WorkflowMetricScan', 'Metric scan mismatch.');
    end
end
finite = isfinite(metric.interval_frame_index) & ...
    isfinite(metric.interval_time_seconds) & isfinite(metric.PQ) & isfinite(metric.PS);
metric = metric(finite, :);
[~, order] = sort(metric.interval_frame_index);
metric = metric(order, :);
end

function [trialTable, scanTable, qc] = build_trial_time_tables(...
        decoded, metric, sessionInfo, cfg, subjectId, scanId)
nBins = cfg.n_time_bins;
dt = cfg.frame_interval_seconds;
pre = cfg.trial_pre_seconds;
post = cfg.trial_post_seconds;
binStarts = -pre + (0:nBins-1) .* dt;
binCenters = binStarts + dt ./ 2;
nTrials = numel(decoded.trials);
nRows = nTrials * nBins;

subjectColumn = repmat({subjectId}, nRows, 1);
scanColumn = repmat({scanId}, nRows, 1);
trialIndex = NaN(nRows, 1);
blockIndex = NaN(nRows, 1);
blockNumber = NaN(nRows, 1);
conditionCode = NaN(nRows, 1);
condition = repmat({''}, nRows, 1);
isPseudo = false(nRows, 1);
trialIndexInBlock = NaN(nRows, 1);
timeIndex = NaN(nRows, 1);
timeSeconds = NaN(nRows, 1);
eventTime = NaN(nRows, 1);
windowStart = NaN(nRows, 1);
windowEnd = NaN(nRows, 1);
validTrial = false(nRows, 1);
coverageFraction = NaN(nRows, 1);
pq = NaN(nRows, 1);
ps = NaN(nRows, 1);

tol = max(1e-9, dt * 1e-4);
binTolerance = 0.5 * dt + tol;
for iTrial = 1:nTrials
    tr = decoded.trials(iTrial);
    rows = (iTrial-1)*nBins + (1:nBins);
    trialIndex(rows) = tr.trial_index;
    blockIndex(rows) = tr.block_index;
    blockNumber(rows) = tr.block_number;
    conditionCode(rows) = tr.condition_code;
    condition(rows) = repmat({condition_name(tr.condition_code)}, nBins, 1);
    isPseudo(rows) = tr.is_pseudo_fixation;
    trialIndexInBlock(rows) = tr.trial_index_in_block;
    timeIndex(rows) = (1:nBins)';
    timeSeconds(rows) = binCenters(:);

    eventTimeThis = (tr.event_sample - sessionInfo.session_start_sample_raw) ./ ...
        sessionInfo.recording_fsample_hz;
    windowStartThis = (tr.window_start_sample - sessionInfo.session_start_sample_raw) ./ ...
        sessionInfo.recording_fsample_hz;
    windowEndThis = (tr.window_end_sample - sessionInfo.session_start_sample_raw) ./ ...
        sessionInfo.recording_fsample_hz;
    parentStart = (tr.parent_start_sample - sessionInfo.session_start_sample_raw) ./ ...
        sessionInfo.recording_fsample_hz;
    parentEnd = (tr.parent_end_sample - sessionInfo.session_start_sample_raw) ./ ...
        sessionInfo.recording_fsample_hz;
    eventTime(rows) = eventTimeThis;
    windowStart(rows) = windowStartThis;
    windowEnd(rows) = windowEndThis;

    if cfg.require_parent_block_containment
        selectStart = max(windowStartThis, parentStart);
        selectEnd = min(windowEndThis, parentEnd);
    else
        selectStart = windowStartThis;
        selectEnd = windowEndThis;
    end
    pqThis = NaN(nBins, 1);
    psThis = NaN(nBins, 1);
    if selectEnd > selectStart
        in = metric.interval_time_seconds >= selectStart - tol & ...
            metric.interval_time_seconds + dt <= selectEnd + tol;
        selected = find(in);
        if ~isempty(selected)
            relativeStart = metric.interval_time_seconds(selected) - eventTimeThis;
            bin = round((relativeStart + pre) ./ dt) + 1;
            for iSelected = 1:numel(selected)
                iBin = bin(iSelected);
                if iBin < 1 || iBin > nBins
                    continue
                end
                % Raw trigger samples need not fall exactly on the 20-ms
                % derivative grid.  Assign to the nearest bin and require the
                % offset to be within half a bin; this preserves alignment
                % without discarding every event because of sample rounding.
                if abs(relativeStart(iSelected) - binStarts(iBin)) > binTolerance
                    continue
                end
                % The first interval wins if a malformed input repeats a bin;
                % normal P04-EXP-003 tables have one interval per bin.
                if ~isfinite(pqThis(iBin))
                    pqThis(iBin) = metric.PQ(selected(iSelected));
                    psThis(iBin) = metric.PS(selected(iSelected));
                end
            end
        end
    end
    present = isfinite(pqThis) & isfinite(psThis);
    coverage = nnz(present) * dt ./ (pre + post);
    isValid = coverage >= cfg.min_trial_coverage_fraction && any(present);
    if ~isValid
        pqThis(:) = NaN;
        psThis(:) = NaN;
    end
    validTrial(rows) = isValid;
    coverageFraction(rows) = coverage;
    pq(rows) = pqThis;
    ps(rows) = psThis;
end

trialTable = table(subjectColumn, scanColumn, trialIndex, blockIndex, ...
    blockNumber, conditionCode, condition, isPseudo, trialIndexInBlock, ...
    timeIndex, timeSeconds, eventTime, windowStart, windowEnd, validTrial, ...
    coverageFraction, pq, ps, 'VariableNames', ...
    {'subject_id','scan_id','trial_index','block_index','source_block_number', ...
    'condition_code','condition','is_pseudo_fixation','trial_index_in_block', ...
    'time_index','time_seconds','event_time_seconds','window_start_seconds', ...
    'window_end_seconds','valid_trial','trial_coverage_fraction','PQ','PS'});

scanTable = aggregate_scan_time(trialTable, cfg, subjectId, scanId, binCenters);
qc = trial_time_qc(trialTable, decoded, cfg);
end

function scanTable = aggregate_scan_time(trialTable, cfg, subjectId, scanId, binCenters)
nBins = cfg.n_time_bins;
nRows = 3 * nBins;
subjectColumn = repmat({subjectId}, nRows, 1);
scanColumn = repmat({scanId}, nRows, 1);
conditionCode = NaN(nRows, 1);
condition = repmat({''}, nRows, 1);
timeIndex = NaN(nRows, 1);
timeSeconds = NaN(nRows, 1);
nValidTrials = zeros(nRows, 1);
nValidPQ = zeros(nRows, 1);
nValidPS = zeros(nRows, 1);
pqMean = NaN(nRows, 1);
psMean = NaN(nRows, 1);
validPQ = false(nRows, 1);
validPS = false(nRows, 1);

row = 0;
for code = 1:3
    for t = 1:nBins
        row = row + 1;
        conditionCode(row) = code;
        condition(row) = {condition_name(code)};
        timeIndex(row) = t;
        timeSeconds(row) = binCenters(t);
        base = trialTable.condition_code == code & ...
            trialTable.time_index == t & trialTable.valid_trial;
        nValidTrials(row) = nnz(base);
        ixPQ = base & isfinite(trialTable.PQ);
        ixPS = base & isfinite(trialTable.PS);
        nValidPQ(row) = nnz(ixPQ);
        nValidPS(row) = nnz(ixPS);
        if nValidPQ(row) >= cfg.min_trials_per_timepoint
            pqMean(row) = mean(trialTable.PQ(ixPQ));
            validPQ(row) = true;
        end
        if nValidPS(row) >= cfg.min_trials_per_timepoint
            psMean(row) = mean(trialTable.PS(ixPS));
            validPS(row) = true;
        end
    end
end
scanTable = table(subjectColumn, scanColumn, conditionCode, condition, ...
    timeIndex, timeSeconds, nValidTrials, nValidPQ, nValidPS, pqMean, psMean, ...
    validPQ, validPS, 'VariableNames', ...
    {'subject_id','scan_id','condition_code','condition','time_index', ...
    'time_seconds','n_valid_trials','n_valid_PQ','n_valid_PS','PQ_mean', ...
    'PS_mean','valid_PQ','valid_PS'});
end

function qc = trial_time_qc(trialTable, decoded, cfg)
trialIds = unique(trialTable.trial_index);
validByTrial = false(numel(trialIds), 1);
coverageByTrial = NaN(numel(trialIds), 1);
for i = 1:numel(trialIds)
    ix = trialTable.trial_index == trialIds(i);
    validByTrial(i) = any(trialTable.valid_trial(ix));
    coverageByTrial(i) = max(trialTable.trial_coverage_fraction(ix));
end
qc = struct;
qc.n_trials_total = numel(decoded.trials);
qc.n_valid_trials = nnz(validByTrial);
qc.n_trial_time_rows = height(trialTable);
qc.n_time_bins = cfg.n_time_bins;
qc.time_window_seconds = [-cfg.trial_pre_seconds, cfg.trial_post_seconds];
qc.coverage_threshold = cfg.min_trial_coverage_fraction;
qc.coverage_median = safe_median(coverageByTrial);
qc.coverage_min = safe_min(coverageByTrial);
qc.condition_trial_counts = struct(...
    'zero_back', nnz(trialTable.condition_code == 1 & ...
        trialTable.time_index == 1), ...
    'two_back', nnz(trialTable.condition_code == 2 & ...
        trialTable.time_index == 1), ...
    'fixation', nnz(trialTable.condition_code == 3 & ...
        trialTable.time_index == 1));
qc.condition_valid_trials = struct(...
    'zero_back', nnz(trialTable.condition_code == 1 & ...
        trialTable.time_index == 1 & trialTable.valid_trial), ...
    'two_back', nnz(trialTable.condition_code == 2 & ...
        trialTable.time_index == 1 & trialTable.valid_trial), ...
    'fixation', nnz(trialTable.condition_code == 3 & ...
        trialTable.time_index == 1 & trialTable.valid_trial));
end

function validate_scan_bounds(decoded, sessionInfo, subjectId, scanId)
if decoded.scan_start_sample ~= double(sessionInfo.session_start_sample_raw) || ...
        decoded.scan_end_sample ~= double(sessionInfo.session_end_sample_raw)
    error('MEG:WorkflowBounds', ...
        'Event/session bounds disagree for %s scan %s.', subjectId, scanId);
end
end

function rawFile = resolve_raw_file(cfg, subjectId, scanId)
rawFile = char(cfg.raw_datafile_template);
rawFile = strrep(rawFile, '{raw_root}', char(cfg.raw_root));
rawFile = strrep(rawFile, '{subject_id}', subjectId);
rawFile = strrep(rawFile, '{scan_id}', scanId);
if isfield(cfg, 'workflow')
    rawFile = strrep(rawFile, '{workflow}', char(cfg.workflow));
end
end

function d = resolve_session_directory(cfg, subjectId, scanId)
name = sprintf('%s_MEG_%s-Wrkmem_task_scan_sessionlcmv_fixedori_power_dt20ms', ...
    subjectId, scanId);
d = fullfile(cfg.processed_root, subjectId, 'Wrkmem', ...
    'srcsessionlcmv_power_exp003', name);
end

function scanIds = normalize_scan_ids(scanIds)
if isnumeric(scanIds)
    scanIds = arrayfun(@num2str, scanIds(:)', 'UniformOutput', false);
elseif ischar(scanIds)
    scanIds = {scanIds};
elseif isstring(scanIds)
    scanIds = cellstr(scanIds(:)');
elseif ~iscell(scanIds)
    error('MEG:WorkflowScanType', 'Invalid scan_ids type.');
end
for i = 1:numel(scanIds)
    scanIds{i} = validate_identifier(scanIds{i}, 'scanId');
end
end

function identifier = validate_identifier(identifier, name)
identifier = char(string(identifier));
if isempty(regexp(identifier, '^[A-Za-z0-9_-]+$', 'once'))
    error('MEG:Identifier', '%s contains unsupported characters.', name);
end
end

function root = resolve_root(cfg, envField, valueField)
root = getenv(char(cfg.(envField)));
if isempty(root), root = char(cfg.(valueField)); end
end

function name = condition_name(code)
if code == 1
    name = '0-back';
elseif code == 2
    name = '2-back';
elseif code == 3
    name = 'fixation';
else
    name = 'unknown';
end
end

function value = safe_median(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = median(x); end
end

function value = safe_min(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = min(x); end
end

function ids = identifiers_to_cellstr(x)
if isnumeric(x)
    ids = arrayfun(@(v) num2str(v, '%.0f'), x(:), 'UniformOutput', false);
elseif iscell(x)
    ids = cellfun(@char, x(:), 'UniformOutput', false);
else
    ids = cellstr(x);
end
end

function write_json(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('MEG:ManifestWrite', 'Cannot open %s.', filename); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value), 'char');
clear cleanup
end
