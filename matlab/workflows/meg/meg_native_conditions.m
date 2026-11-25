function summary = meg_native_conditions(configFile, subjectId, scanIds)
%MEG_NATIVE_CONDITIONS Align P04-EXP-003 PQ/PS metrics to HCP WM events.
%
% One array task processes one participant.  Raw trigger events are decoded
% with the official HCP trial function; no condition labels are inferred from
% the averaged or velocity derivatives.

if nargin < 2 || nargin > 3
    error('MEG:WorkflowArguments', ...
        'Usage: meg_native_conditions(configFile, subjectId, [scanIds])');
end
cfg = load_json_config(configFile);
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
    'raw_file', '', 'metric_file', '', 'block_file', '', 'trial_file', '', ...
    'n_blocks', 0, 'n_trials', 0, 'n_valid_blocks', 0, ...
    'n_valid_trials', 0, 'error_identifier', '', 'error_message', ''), ...
    1, numel(scanIds));

for iScan = 1:numel(scanIds)
    scanId = scanIds{iScan};
    scanSummary(iScan).scan_id = scanId;
    scanResultRoot = fullfile(subjectResultRoot, ['scan-' scanId]);
    if ~isfolder(scanResultRoot), mkdir(scanResultRoot); end
    try
        rawFile = resolve_raw_file(cfg, subjectId, scanId);
        metricDir = resolve_metric_directory(cfg, subjectId, scanId);
        metricFile = fullfile(metricDir, 'pq_ps_per_interval.csv');
        sessionDir = resolve_session_directory(cfg, subjectId, scanId);
        metadataFile = fullfile(sessionDir, 'session_metadata.mat');
        required = {rawFile, metricFile, metadataFile, ...
            fullfile(sessionDir, 'session_manifest.json')};
        for iFile = 1:numel(required)
            if ~isfile(required{iFile})
                error('MEG:WorkflowInputMissing', 'Required input not found: %s', required{iFile});
            end
        end

        loaded = load(metadataFile, 'metadata');
        if ~isfield(loaded, 'metadata') || ~isfield(loaded.metadata, 'session_input_info')
            error('MEG:WorkflowMetadata', 'session_metadata.mat lacks session_input_info.');
        end
        sessionInfo = loaded.metadata.session_input_info;
        metric = read_metric_table(metricFile, subjectId, scanId);
        decoded = meg_decode_events(rawFile, cfg.megconnectome_root, ...
            cfg.trial_pre_seconds, cfg.trial_post_seconds);
        validate_scan_bounds(decoded, sessionInfo, subjectId, scanId);

        [blockTable, blockQC] = aggregate_blocks(decoded, metric, sessionInfo, cfg, subjectId, scanId);
        [trialTable, trialQC] = aggregate_trials(decoded, metric, sessionInfo, cfg, subjectId, scanId);
        writetable(blockTable, fullfile(scanResultRoot, 'block_metrics.csv'));
        writetable(trialTable, fullfile(scanResultRoot, 'trial_metrics.csv'));

        alignment = struct;
        alignment.schema_version = 1;
        alignment.status = 'complete';
        alignment.workflow = 'meg_native_conditions';
        alignment.parent_metric_experiment = char(cfg.metric_results_subdir);
        alignment.subject_id = subjectId;
        alignment.scan_id = scanId;
        alignment.raw_file = rawFile;
        alignment.metric_file = metricFile;
        alignment.session_metadata_file = metadataFile;
        alignment.recording_fsample_hz = decoded.recording_fsample_hz;
        alignment.recording_samples = decoded.recording_samples;
        alignment.scan_start_sample = decoded.scan_start_sample;
        alignment.scan_end_sample = decoded.scan_end_sample;
        alignment.n_official_blocks = numel(decoded.blocks);
        alignment.n_official_trials = numel(decoded.trials);
        alignment.n_block_rows = height(blockTable);
        alignment.n_trial_rows = height(trialTable);
        alignment.n_valid_blocks = nnz(blockTable.valid_window);
        alignment.n_valid_trials = nnz(trialTable.valid_window);
        alignment.block_qc = blockQC;
        alignment.trial_qc = trialQC;
        alignment.decoder_warning_count = decoded.warning_count;
        alignment.created = datestr(now, 30);
        write_json(fullfile(scanResultRoot, 'alignment.json'), alignment);
        save(fullfile(scanResultRoot, 'alignment.mat'), 'alignment', 'blockQC', 'trialQC');

        scanSummary(iScan).status = 'complete';
        scanSummary(iScan).raw_file = rawFile;
        scanSummary(iScan).metric_file = metricFile;
        scanSummary(iScan).block_file = fullfile(scanResultRoot, 'block_metrics.csv');
        scanSummary(iScan).trial_file = fullfile(scanResultRoot, 'trial_metrics.csv');
        scanSummary(iScan).n_blocks = height(blockTable);
        scanSummary(iScan).n_trials = height(trialTable);
        scanSummary(iScan).n_valid_blocks = nnz(blockTable.valid_window);
        scanSummary(iScan).n_valid_trials = nnz(trialTable.valid_window);
    catch ME
        scanSummary(iScan).status = 'failed';
        scanSummary(iScan).error_identifier = ME.identifier;
        scanSummary(iScan).error_message = getReport(ME, 'extended', 'hyperlinks', 'off');
        failure = struct('status', 'failed', 'subject_id', subjectId, ...
            'scan_id', scanId, 'error_identifier', ME.identifier, ...
            'error_message', scanSummary(iScan).error_message, ...
            'created', datestr(now, 30));
        write_json(fullfile(scanResultRoot, 'alignment.json'), failure);
    end
end

summary = struct;
summary.schema_version = 1;
summary.workflow = 'meg_native_conditions';
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
    error('MEG:WorkflowSubjectFailed', 'meg_native_conditions failed for subject %s.', subjectId);
end
end

function cfg = load_json_config(configFile)
if nargin ~= 1 || ~isfile(configFile)
    error('MEG:WorkflowConfig', 'Configuration file not found: %s', char(configFile));
end
cfg = jsondecode(fileread(configFile));
required = {'workflow','raw_root_env','raw_root','processed_root_env', ...
    'processed_root','results_root_env','results_root','fieldtrip_root_env', ...
    'fieldtrip_root','megconnectome_root_env','megconnectome_root', ...
    'raw_datafile_template','metric_results_subdir','results_subdir', ...
    'scan_ids','frame_interval_seconds','trial_pre_seconds','trial_post_seconds', ...
    'min_window_coverage_fraction','require_parent_block_containment'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i}), error('MEG:WorkflowConfigField', 'Missing field %s.', required{i}); end
end
cfg.raw_root = resolve_root(cfg, 'raw_root_env', 'raw_root');
cfg.processed_root = resolve_root(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_root(cfg, 'results_root_env', 'results_root');
cfg.fieldtrip_root = resolve_root(cfg, 'fieldtrip_root_env', 'fieldtrip_root');
cfg.megconnectome_root = resolve_root(cfg, 'megconnectome_root_env', 'megconnectome_root');
for f = {'raw_root','processed_root','results_root','fieldtrip_root','megconnectome_root'}
    if ~isfolder(cfg.(f{1})), error('MEG:WorkflowRoot', 'Root not found: %s', cfg.(f{1})); end
end
if cfg.frame_interval_seconds <= 0 || cfg.trial_pre_seconds < 0 || cfg.trial_post_seconds <= 0
    error('MEG:WorkflowConfigValue', 'Invalid time parameters.');
end
end

function metric = read_metric_table(filename, subjectId, scanId)
metric = readtable(filename);
required = {'interval_frame_index','interval_time_seconds','PQ','PS','PQ_curl','PS_curl'};
for i = 1:numel(required)
    if ~ismember(required{i}, metric.Properties.VariableNames)
        error('MEG:WorkflowMetricColumns', 'Missing metric column %s in %s.', required{i}, filename);
    end
end
if ismember('subject_id', metric.Properties.VariableNames)
    ids = identifiers_to_cellstr(metric.subject_id);
    if any(~strcmp(ids, subjectId)), error('MEG:WorkflowMetricSubject', 'Metric subject mismatch.'); end
end
if ismember('scan_id', metric.Properties.VariableNames)
    ids = identifiers_to_cellstr(metric.scan_id);
    if any(~strcmp(ids, scanId)), error('MEG:WorkflowMetricScan', 'Metric scan mismatch.'); end
end
finite = isfinite(metric.interval_frame_index) & isfinite(metric.interval_time_seconds) & ...
    isfinite(metric.PQ) & isfinite(metric.PS) & isfinite(metric.PQ_curl) & isfinite(metric.PS_curl);
metric = metric(finite, :);
[~, order] = sort(metric.interval_frame_index);
metric = metric(order, :);
end

function [blockTable, qc] = aggregate_blocks(decoded, metric, sessionInfo, cfg, subjectId, scanId)
n = numel(decoded.blocks);
rows = repmat(empty_block_row(), n, 1);
for i = 1:n
    b = decoded.blocks(i);
    tStart = (b.start_sample - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    tEnd = (b.end_sample - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    [vals, coverage, in] = select_intervals(metric, tStart, tEnd, cfg.frame_interval_seconds, tEnd - tStart);
    rows(i) = summarize_window(vals, coverage, tStart, tEnd, in, ...
        subjectId, scanId, b.block_number, b.block_index, b.condition_code, ...
        b.image_type, NaN, tStart, tEnd, NaN, false, cfg.min_window_coverage_fraction);
    rows(i).trial_index = NaN;
end
blockTable = block_rows_to_table(rows);
qc = qc_summary(blockTable, 'block');
end

function [trialTable, qc] = aggregate_trials(decoded, metric, sessionInfo, cfg, subjectId, scanId)
n = numel(decoded.trials);
rows = repmat(empty_trial_row(), n, 1);
for i = 1:n
    tr = decoded.trials(i);
    rawStart = tr.window_start_sample;
    rawEnd = tr.window_end_sample;
    windowStart = (rawStart - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    windowEnd = (rawEnd - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    parentStart = (tr.parent_start_sample - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    parentEnd = (tr.parent_end_sample - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    if cfg.require_parent_block_containment
        selectStart = max(windowStart, parentStart);
        selectEnd = min(windowEnd, parentEnd);
    else
        selectStart = windowStart;
        selectEnd = windowEnd;
    end
    [vals, coverage, in] = select_intervals(metric, selectStart, selectEnd, ...
        cfg.frame_interval_seconds, windowEnd - windowStart);
    eventTime = (tr.event_sample - sessionInfo.session_start_sample_raw) / sessionInfo.recording_fsample_hz;
    rows(i) = summarize_window(vals, coverage, selectStart, selectEnd, in, ...
        subjectId, scanId, tr.trial_index, tr.block_index, tr.condition_code, ...
        tr.image_type, tr.trial_index_in_block, parentStart, parentEnd, ...
        NaN, tr.is_pseudo_fixation, cfg.min_window_coverage_fraction);
    rows(i).source_block_number = tr.block_number;
    rows(i).trial_index = tr.trial_index;
    rows(i).window_start_seconds = windowStart;
    rows(i).window_end_seconds = windowEnd;
    rows(i).event_time_seconds = eventTime;
end
trialTable = trial_rows_to_table(rows);
qc = qc_summary(trialTable, 'trial');
end

function [vals, coverage, in] = select_intervals(metric, startTime, endTime, dt, denominator)
if endTime <= startTime
    in = false(height(metric), 1);
else
    % An interval is included only when the complete 20-ms interval lies in
    % the requested window.  The small tolerance absorbs raw-sample rounding.
    tol = max(1e-9, dt * 1e-4);
    in = metric.interval_time_seconds >= startTime - tol & ...
        metric.interval_time_seconds + dt <= endTime + tol;
end
vals = [metric.PQ(in), metric.PS(in), metric.PQ_curl(in), metric.PS_curl(in)];
if denominator <= 0, coverage = 0; else, coverage = size(vals, 1) * dt / denominator; end
coverage = min(1, max(0, coverage));
end

function row = empty_block_row()
row = struct('subject_id', '', 'scan_id', '', 'block_index', NaN, ...
    'source_block_number', NaN, 'condition_code', NaN, 'condition', '', ...
    'image_type', NaN, 'trial_index', NaN, 'window_start_seconds', NaN, ...
    'window_end_seconds', NaN, 'event_time_seconds', NaN, ...
    'parent_start_seconds', NaN, 'parent_end_seconds', NaN, ...
    'n_intervals', 0, 'coverage_fraction', 0, 'valid_window', false, ...
    'PQ_mean', NaN, 'PQ_median', NaN, 'PQ_sd', NaN, ...
    'PS_mean', NaN, 'PS_median', NaN, 'PS_sd', NaN, ...
    'PQ_curl_mean', NaN, 'PS_curl_mean', NaN);
row.trial_index_in_block = NaN;
row.is_pseudo_fixation = false;
end

function row = empty_trial_row()
row = empty_block_row();
end

function row = summarize_window(vals, coverage, startTime, endTime, in, subjectId, scanId, ...
        rowIndex, blockIndex, conditionCode, imageType, trialIndexInBlock, ...
        parentStart, parentEnd, eventTime, isPseudo, minCoverage)
row = empty_block_row();
row.subject_id = subjectId;
row.scan_id = scanId;
row.block_index = blockIndex;
row.source_block_number = rowIndex;
row.condition_code = conditionCode;
row.condition = condition_name(conditionCode);
row.image_type = imageType;
row.trial_index = rowIndex;
row.window_start_seconds = startTime;
row.window_end_seconds = endTime;
row.event_time_seconds = eventTime;
row.parent_start_seconds = parentStart;
row.parent_end_seconds = parentEnd;
row.n_intervals = nnz(in);
row.coverage_fraction = coverage;
row.valid_window = coverage >= minCoverage && row.n_intervals > 0;
if row.n_intervals > 0
    row.PQ_mean = mean(vals(:,1)); row.PQ_median = median(vals(:,1));
    row.PQ_sd = safe_sd(vals(:,1));
    row.PS_mean = mean(vals(:,2)); row.PS_median = median(vals(:,2));
    row.PS_sd = safe_sd(vals(:,2));
    row.PQ_curl_mean = mean(vals(:,3)); row.PS_curl_mean = mean(vals(:,4));
end
if isfield(row, 'trial_index_in_block'), row.trial_index_in_block = trialIndexInBlock; end
if isfield(row, 'is_pseudo_fixation'), row.is_pseudo_fixation = isPseudo; end
end

function tbl = block_rows_to_table(rows)
tbl = struct2table(rows);
% The block index is the temporal index after task and Fix15sec boundaries
% were merged; source_block_number is the official task block number.
end

function tbl = trial_rows_to_table(rows)
tbl = struct2table(rows);
end

function name = condition_name(code)
if code == 1, name = '0-back';
elseif code == 2, name = '2-back';
elseif code == 3, name = 'fixation';
else, name = 'unknown';
end
end

function qc = qc_summary(tbl, kind)
qc = struct;
qc.kind = kind;
qc.n_rows = height(tbl);
qc.n_valid = nnz(tbl.valid_window);
qc.coverage_median = safe_median(tbl.coverage_fraction);
qc.coverage_min = safe_min(tbl.coverage_fraction);
qc.condition_counts = struct('zero_back', nnz(tbl.condition_code == 1), ...
    'two_back', nnz(tbl.condition_code == 2), 'fixation', nnz(tbl.condition_code == 3));
qc.valid_condition_counts = struct('zero_back', nnz(tbl.valid_window & tbl.condition_code == 1), ...
    'two_back', nnz(tbl.valid_window & tbl.condition_code == 2), ...
    'fixation', nnz(tbl.valid_window & tbl.condition_code == 3));
end

function validate_scan_bounds(decoded, sessionInfo, subjectId, scanId)
if decoded.scan_start_sample ~= double(sessionInfo.session_start_sample_raw) || ...
        decoded.scan_end_sample ~= double(sessionInfo.session_end_sample_raw)
    error('MEG:WorkflowBounds', 'Event/session bounds disagree for %s scan %s.', subjectId, scanId);
end
end

function rawFile = resolve_raw_file(cfg, subjectId, scanId)
rawFile = char(cfg.raw_datafile_template);
rawFile = strrep(rawFile, '{raw_root}', char(cfg.raw_root));
rawFile = strrep(rawFile, '{subject_id}', subjectId);
rawFile = strrep(rawFile, '{scan_id}', scanId);
end

function d = resolve_metric_directory(cfg, subjectId, scanId)
scope = 'task_scan';
name = sprintf('%s_MEG_%s-Wrkmem_%s_sessionlcmv_fixedori_power_dt20ms', ...
    subjectId, scanId, scope);
d = fullfile(cfg.results_root, char(cfg.metric_results_subdir), subjectId, ['scan-' scanId]);
% The compact table is in the P04-EXP-003 results directory, not in the
% processed source directory.  Keep this function explicit for auditability.
if ~isfolder(d), error('MEG:WorkflowMetricDir', 'Metric result directory not found: %s', d); end
end

function d = resolve_session_directory(cfg, subjectId, scanId)
name = sprintf('%s_MEG_%s-Wrkmem_task_scan_sessionlcmv_fixedori_power_dt20ms', ...
    subjectId, scanId);
d = fullfile(cfg.processed_root, subjectId, 'Wrkmem', ...
    'srcsessionlcmv_power_exp003', name);
end

function scanIds = normalize_scan_ids(scanIds)
if isnumeric(scanIds), scanIds = arrayfun(@num2str, scanIds(:)', 'UniformOutput', false);
elseif ischar(scanIds), scanIds = {scanIds};
elseif isstring(scanIds), scanIds = cellstr(scanIds(:)');
elseif ~iscell(scanIds), error('MEG:WorkflowScanType', 'Invalid scan_ids type.');
end
for i = 1:numel(scanIds), scanIds{i} = validate_identifier(scanIds{i}, 'scanId'); end
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

function value = safe_sd(x)
x = x(isfinite(x));
if numel(x) < 2, value = NaN; else, value = std(x); end
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
