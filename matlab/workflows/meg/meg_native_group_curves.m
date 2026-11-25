function result = meg_native_group_curves(configFile)
% Participant-level means and SEM for native MEG trial curves.
% Equal-weight scan aggregation is preserved; no timepoint inference.

if nargin ~= 1
    error('MEG:WorkflowGroupArguments', ...
        'Usage: meg_native_group_curves(configFile)');
end
cfg = load_config(configFile);
projectRoot = fileparts(fileparts(char(configFile)));
manifestPath = meg_subject_manifest(cfg, projectRoot);
subjects = read_subject_manifest(manifestPath);
if isempty(subjects)
    error('MEG:WorkflowGroupSubjects', 'No subjects found in %s.', manifestPath);
end
scanIds = normalize_scan_ids(cfg.scan_ids);
nSubjects = numel(subjects);
nScans = numel(scanIds);
nBins = cfg.n_time_bins;
timeSeconds = -cfg.trial_pre_seconds + ((1:nBins) - 0.5) .* ...
    cfg.frame_interval_seconds;

resultsRoot = resolve_root(cfg, 'results_root_env', 'results_root');
outputRoot = fullfile(resultsRoot, char(cfg.results_subdir), 'group_stats');
if ~isfolder(outputRoot), mkdir(outputRoot); end

scanPQ = NaN(nSubjects, nScans, 3, nBins);
scanPS = NaN(nSubjects, nScans, 3, nBins);
scanTrialN = zeros(nSubjects, nScans, 3, nBins);
failures = struct('subject_id', {}, 'scan_id', {}, 'reason', {});
loadedScans = 0;
for iSubject = 1:nSubjects
    subjectId = subjects{iSubject};
    for iScan = 1:nScans
        scanId = scanIds{iScan};
        scanFile = fullfile(resultsRoot, char(cfg.results_subdir), ...
            subjectId, ['scan-' scanId], 'scan_time_means.csv');
        if ~isfile(scanFile)
            failures(end+1) = failure_row(subjectId, scanId, ...
                'missing scan_time_means.csv'); %#ok<AGROW>
            continue
        end
        try
            tbl = readtable(scanFile);
            validate_scan_table(tbl, subjectId, scanId, cfg);
            for code = 1:3
                for t = 1:nBins
                    ix = tbl.condition_code == code & tbl.time_index == t;
                    if nnz(ix) ~= 1
                        error('MEG:WorkflowScanRows', ...
                            'Expected one row for condition %d/time %d.', code, t);
                    end
                    row = find(ix, 1, 'first');
                    scanTrialN(iSubject, iScan, code, t) = tbl.n_valid_trials(row);
                    if tbl.valid_PQ(row) && isfinite(tbl.PQ_mean(row))
                        scanPQ(iSubject, iScan, code, t) = tbl.PQ_mean(row);
                    end
                    if tbl.valid_PS(row) && isfinite(tbl.PS_mean(row))
                        scanPS(iSubject, iScan, code, t) = tbl.PS_mean(row);
                    end
                end
            end
            loadedScans = loadedScans + 1;
        catch ME
            failures(end+1) = failure_row(subjectId, scanId, ME.message); %#ok<AGROW>
        end
    end
end

audit = struct;
audit.schema_version = 1;
audit.workflow = 'meg_native_trial_curves';
audit.status = 'running';
audit.subject_manifest = manifestPath;
audit.expected_subjects = nSubjects;
audit.expected_scans = nScans;
audit.loaded_scan_tables = loadedScans;
audit.failure_count = numel(failures);
audit.failures = failures;
audit.trial_pre_seconds = cfg.trial_pre_seconds;
audit.trial_post_seconds = cfg.trial_post_seconds;
audit.frame_interval_seconds = cfg.frame_interval_seconds;
audit.n_time_bins = nBins;
audit.created = datestr(now, 30);
write_json(fullfile(outputRoot, 'trial_time_manifest.json'), audit);
if numel(failures) > 0 || loadedScans ~= nSubjects * nScans
    audit.status = 'incomplete';
    write_json(fullfile(outputRoot, 'trial_time_manifest.json'), audit);
    error('MEG:WorkflowGroupIncomplete', ...
        'Subject-level inputs incomplete (%d failures; loaded %d of %d scans).', ...
        numel(failures), loadedScans, nSubjects*nScans);
end

% Equal-weight scan aggregation within participant.  The complete-scan flags
% are retained in the participant table so the inferential denominator is
% auditable at each time point.
participantPQ = NaN(nSubjects, 3, nBins);
participantPS = NaN(nSubjects, 3, nBins);
nScanPQ = zeros(nSubjects, 3, nBins);
nScanPS = zeros(nSubjects, 3, nBins);
for iSubject = 1:nSubjects
    for code = 1:3
        for t = 1:nBins
            values = squeeze(scanPQ(iSubject, :, code, t));
            finite = isfinite(values);
            nScanPQ(iSubject, code, t) = nnz(finite);
            if any(finite), participantPQ(iSubject, code, t) = mean(values(finite)); end
            values = squeeze(scanPS(iSubject, :, code, t));
            finite = isfinite(values);
            nScanPS(iSubject, code, t) = nnz(finite);
            if any(finite), participantPS(iSubject, code, t) = mean(values(finite)); end
        end
    end
end

participantTable = make_participant_table(subjects, timeSeconds, ...
    participantPQ, participantPS, nScanPQ, nScanPS, nScans);
writetable(participantTable, fullfile(outputRoot, ...
    'trial_time_participant_means.csv'));

fixationReferencePQ = finite_mean(participantPQ(:, 3, :));
fixationReferencePS = finite_mean(participantPS(:, 3, :));
fixationReferencePQSem = finite_sem(participantPQ(:, 3, :));
fixationReferencePSSem = finite_sem(participantPS(:, 3, :));
groupTable = make_group_curve_table(timeSeconds, participantPQ, participantPS, ...
    fixationReferencePQ, fixationReferencePS, fixationReferencePQSem, ...
    fixationReferencePSSem);
writetable(groupTable, fullfile(outputRoot, 'trial_time_group_curves.csv'));

result = struct;
result.schema_version = 1;
result.workflow = 'meg_native_trial_curves';
result.status = 'complete';
result.subject_count = nSubjects;
result.scan_count = nScans;
result.time_window_seconds = [-cfg.trial_pre_seconds, cfg.trial_post_seconds];
result.frame_interval_seconds = cfg.frame_interval_seconds;
result.n_time_bins = nBins;
result.fixation_reference_PQ = fixationReferencePQ;
result.fixation_reference_PS = fixationReferencePS;
result.created = datestr(now, 30);
save(fullfile(outputRoot, 'trial_time_group.mat'), 'result', ...
    'participantPQ', 'participantPS', 'nScanPQ', 'nScanPS', ...
    'groupTable', '-v7.3');

audit.status = 'complete';
audit.output_root = outputRoot;
audit.participant_rows = height(participantTable);
audit.group_curve_rows = height(groupTable);
audit.created = result.created;
write_json(fullfile(outputRoot, 'trial_time_manifest.json'), audit);
end

function cfg = load_config(configFile)
cfg = jsondecode(fileread(configFile));
required = {'raw_root_env','raw_root','processed_root_env','processed_root', ...
    'results_root_env','results_root','fieldtrip_root_env','fieldtrip_root', ...
    'megconnectome_root_env','megconnectome_root','results_subdir', ...
    'subject_manifest_relative','scan_ids','frame_interval_seconds', ...
    'n_time_bins','trial_pre_seconds','trial_post_seconds', ...
    'min_trials_per_timepoint'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i})
        error('MEG:WorkflowGroupConfig', 'Missing field %s.', required{i});
    end
end
cfg.raw_root = resolve_root(cfg, 'raw_root_env', 'raw_root');
cfg.processed_root = resolve_root(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_root(cfg, 'results_root_env', 'results_root');
cfg.fieldtrip_root = resolve_root(cfg, 'fieldtrip_root_env', 'fieldtrip_root');
cfg.megconnectome_root = resolve_root(cfg, 'megconnectome_root_env', 'megconnectome_root');
if cfg.frame_interval_seconds <= 0 || cfg.n_time_bins < 1 || ...
        cfg.trial_pre_seconds < 0 || cfg.trial_post_seconds <= 0
    error('MEG:WorkflowGroupConfigValue', 'Invalid time parameters.');
end
expectedBins = round((cfg.trial_pre_seconds + cfg.trial_post_seconds) / ...
    cfg.frame_interval_seconds);
if expectedBins ~= cfg.n_time_bins
    error('MEG:WorkflowGroupConfigBins', 'n_time_bins does not match the time window.');
end
end

function validate_scan_table(tbl, subjectId, scanId, cfg)
required = {'subject_id','scan_id','condition_code','time_index', ...
    'time_seconds','n_valid_trials','valid_PQ','valid_PS','PQ_mean','PS_mean'};
for i = 1:numel(required)
    if ~ismember(required{i}, tbl.Properties.VariableNames)
        error('MEG:WorkflowScanColumns', 'Scan table lacks %s.', required{i});
    end
end
ids = identifiers_to_cellstr(tbl.subject_id);
scans = identifiers_to_cellstr(tbl.scan_id);
if any(~strcmp(ids, subjectId)) || any(~strcmp(scans, scanId))
    error('MEG:WorkflowScanIdentity', 'Scan table identity mismatch.');
end
if height(tbl) ~= 3 * cfg.n_time_bins
    error('MEG:WorkflowScanRowCount', 'Unexpected scan table row count.');
end
if any(~ismember(tbl.condition_code, [1 2 3])) || ...
        any(tbl.time_index < 1) || any(tbl.time_index > cfg.n_time_bins)
    error('MEG:WorkflowScanIndex', 'Unexpected condition/time index.');
end
end

function tbl = make_participant_table(subjects, timeSeconds, participantPQ, ...
        participantPS, nScanPQ, nScanPS, nScans)
nSubjects = numel(subjects);
nBins = numel(timeSeconds);
nRows = nSubjects * 3 * nBins;
subjectColumn = repmat({''}, nRows, 1);
conditionCode = NaN(nRows, 1);
condition = repmat({''}, nRows, 1);
timeIndex = NaN(nRows, 1);
timeColumn = NaN(nRows, 1);
nScansPQ = zeros(nRows, 1);
nScansPS = zeros(nRows, 1);
pq = NaN(nRows, 1);
ps = NaN(nRows, 1);
completePQ = false(nRows, 1);
completePS = false(nRows, 1);
row = 0;
for iSubject = 1:nSubjects
    for code = 1:3
        for t = 1:nBins
            row = row + 1;
            subjectColumn(row) = subjects(iSubject);
            conditionCode(row) = code;
            condition(row) = {condition_name(code)};
            timeIndex(row) = t;
            timeColumn(row) = timeSeconds(t);
            nScansPQ(row) = nScanPQ(iSubject, code, t);
            nScansPS(row) = nScanPS(iSubject, code, t);
            pq(row) = participantPQ(iSubject, code, t);
            ps(row) = participantPS(iSubject, code, t);
            completePQ(row) = nScansPQ(row) == nScans;
            completePS(row) = nScansPS(row) == nScans;
        end
    end
end
tbl = table(subjectColumn, conditionCode, condition, timeIndex, timeColumn, ...
    nScansPQ, nScansPS, pq, ps, completePQ, completePS, ...
    'VariableNames', {'subject_id','condition_code','condition','time_index', ...
    'time_seconds','n_scans_PQ','n_scans_PS','PQ_mean','PS_mean', ...
    'complete_PQ','complete_PS'});
end

function tbl = make_group_curve_table(timeSeconds, participantPQ, participantPS, ...
        fixPQ, fixPS, fixPQSem, fixPSSem)
nBins = numel(timeSeconds);
nRows = 3 * nBins;
conditionCode = NaN(nRows, 1);
condition = repmat({''}, nRows, 1);
timeIndex = NaN(nRows, 1);
timeColumn = NaN(nRows, 1);
nPQ = zeros(nRows, 1); nPS = zeros(nRows, 1);
pqMean = NaN(nRows, 1); pqSem = NaN(nRows, 1);
psMean = NaN(nRows, 1); psSem = NaN(nRows, 1);
fixPQColumn = repmat(fixPQ, nRows, 1);
fixPSColumn = repmat(fixPS, nRows, 1);
fixPQSemColumn = repmat(fixPQSem, nRows, 1);
fixPSSemColumn = repmat(fixPSSem, nRows, 1);
row = 0;
for code = 1:3
    for t = 1:nBins
        row = row + 1;
        conditionCode(row) = code;
        condition(row) = {condition_name(code)};
        timeIndex(row) = t;
        timeColumn(row) = timeSeconds(t);
        values = participantPQ(:, code, t);
        nPQ(row) = nnz(isfinite(values));
        pqMean(row) = finite_mean(values);
        pqSem(row) = finite_sem(values);
        values = participantPS(:, code, t);
        nPS(row) = nnz(isfinite(values));
        psMean(row) = finite_mean(values);
        psSem(row) = finite_sem(values);
    end
end
tbl = table(conditionCode, condition, timeIndex, timeColumn, nPQ, pqMean, ...
    pqSem, nPS, psMean, psSem, fixPQColumn, fixPQSemColumn, fixPSColumn, ...
    fixPSSemColumn, 'VariableNames', {'condition_code','condition','time_index', ...
    'time_seconds','n_participants_PQ','PQ_mean','PQ_sem','n_participants_PS', ...
    'PS_mean','PS_sem','fixation_reference_PQ','fixation_reference_PQ_sem', ...
    'fixation_reference_PS','fixation_reference_PS_sem'});
end

function value = finite_mean(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = mean(x(:)); end
end

function value = finite_median(x)
x = x(isfinite(x));
if isempty(x), value = NaN; else, value = median(x(:)); end
end

function value = finite_sem(x)
x = x(isfinite(x));
if numel(x) < 2, value = NaN; else, value = std(x(:)) ./ sqrt(numel(x)); end
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

function failures = failure_row(subjectId, scanId, reason)
failures = struct('subject_id',subjectId,'scan_id',scanId,'reason',reason);
end

function subjects = read_subject_manifest(filename)
if ~isfile(filename)
    error('MEG:WorkflowSubjectManifest', 'Subject manifest not found: %s', filename);
end
fid = fopen(filename, 'r');
cleanup = onCleanup(@() fclose(fid));
data = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
subjects = data{1};
subjects = subjects(~cellfun(@isempty, subjects));
for i = 1:numel(subjects)
    subjects{i} = strtrim(subjects{i});
    if isempty(regexp(subjects{i}, '^[0-9]+$', 'once'))
        error('MEG:WorkflowSubjectManifest', ...
            'Invalid subject identifier: %s', subjects{i});
    end
end
subjects = unique(subjects, 'stable');
clear cleanup
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

function root = resolve_root(cfg, envField, valueField)
root = getenv(char(cfg.(envField)));
if isempty(root), root = char(cfg.(valueField)); end
end

function write_json(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('MEG:ManifestWrite', 'Cannot open %s.', filename); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value), 'char');
clear cleanup
end
