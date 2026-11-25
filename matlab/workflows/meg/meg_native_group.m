function result = meg_native_group(configFile)
%MEG_NATIVE_GROUP Participant-aware cohort PQ/PS statistics.
%
% Subject/scan tables are produced by meg_native_conditions.  This entry
% point only aggregates those tables and never treats 20-ms intervals or
% repeated trials as independent population observations.

if nargin ~= 1
    error('MEG:WorkflowGroupArguments', 'Usage: meg_native_group(configFile)');
end
cfg = load_config(configFile);
projectRoot = fileparts(fileparts(char(configFile)));
manifestPath = meg_subject_manifest(cfg, projectRoot);
subjects = read_subject_manifest(manifestPath);
if isempty(subjects)
    error('MEG:WorkflowGroupSubjects', 'No subjects found in %s.', manifestPath);
end
scanIds = normalize_scan_ids(cfg.scan_ids);
resultsRoot = cfg.results_root;
outputRoot = fullfile(resultsRoot, char(cfg.results_subdir), 'group_stats');
if ~isfolder(outputRoot), mkdir(outputRoot); end

blockTables = cell(0, 1);
trialTables = cell(0, 1);
failures = struct('subject_id', {}, 'scan_id', {}, 'reason', {});
for iSubject = 1:numel(subjects)
    subjectId = subjects{iSubject};
    for iScan = 1:numel(scanIds)
        scanId = scanIds{iScan};
        scanRoot = fullfile(resultsRoot, char(cfg.results_subdir), subjectId, ['scan-' scanId]);
        blockFile = fullfile(scanRoot, 'block_metrics.csv');
        trialFile = fullfile(scanRoot, 'trial_metrics.csv');
        if ~isfile(blockFile) || ~isfile(trialFile)
            failures(end+1) = failure_row(subjectId, scanId, 'missing subject-level table'); %#ok<AGROW>
            continue
        end
        try
            bt = readtable(blockFile);
            tt = readtable(trialFile);
            validate_subject_table(bt, subjectId, scanId, 'block');
            validate_subject_table(tt, subjectId, scanId, 'trial');
            blockTables{end+1} = bt; %#ok<AGROW>
            trialTables{end+1} = tt; %#ok<AGROW>
        catch ME
            failures(end+1) = failure_row(subjectId, scanId, ME.message); %#ok<AGROW>
        end
    end
end

audit = struct;
audit.schema_version = 1;
audit.workflow = 'meg_native_conditions';
audit.status = 'running';
audit.subject_manifest = manifestPath;
audit.expected_subjects = numel(subjects);
audit.expected_scans = numel(scanIds);
audit.loaded_scan_tables = numel(blockTables);
audit.failure_count = numel(failures);
audit.failures = failures;
audit.created = datestr(now, 30);
write_json(fullfile(outputRoot, 'group_stats_manifest.json'), audit);
if numel(failures) > 0 || numel(blockTables) ~= numel(subjects) * numel(scanIds)
    audit.status = 'incomplete';
    write_json(fullfile(outputRoot, 'group_stats_manifest.json'), audit);
    error('MEG:WorkflowGroupIncomplete', ...
        'Subject-level inputs incomplete (%d failures; loaded %d of %d scans).', ...
        numel(failures), numel(blockTables), numel(subjects) * numel(scanIds));
end

allBlocks = vertcat(blockTables{:});
allTrials = vertcat(trialTables{:});
blockScan = struct2table(aggregate_scan_means(allBlocks, subjects, scanIds));
trialScan = struct2table(aggregate_scan_means(allTrials, subjects, scanIds));
blockParticipant = struct2table(aggregate_participant_means(blockScan, subjects, scanIds));
trialParticipant = struct2table(aggregate_participant_means(trialScan, subjects, scanIds));
writetable(blockScan, fullfile(outputRoot, 'block_scan_means.csv'));
writetable(trialScan, fullfile(outputRoot, 'trial_scan_means.csv'));
writetable(blockParticipant, fullfile(outputRoot, 'block_participant_means.csv'));
writetable(trialParticipant, fullfile(outputRoot, 'trial_participant_means.csv'));

% Retain the six original block comparisons and bootstrap confidence intervals.
% The Python manuscript analysis independently verifies their p values.
blockStats = analyze_block_level(blockParticipant,subjects,cfg);
trialStats = describe_conditions(trialParticipant,subjects,'trial');
write_struct_table(blockStats.pairwise,fullfile(outputRoot,'block_pairwise.csv'));
write_struct_table([blockStats.descriptive; trialStats],fullfile(outputRoot,'condition_descriptives.csv'));

result = struct;
result.schema_version = 1;
result.workflow = 'meg_native_conditions';
result.status = 'complete';
result.subject_count = numel(subjects);
result.scan_count = numel(scanIds);
result.block = blockStats;
result.trial = trialStats;
result.created = datestr(now, 30);
save(fullfile(outputRoot, 'group_stats.mat'), 'result', 'blockScan', ...
    'trialScan', 'blockParticipant', 'trialParticipant', '-v7.3');

audit.status = 'complete';
audit.output_root = outputRoot;
audit.block_rows = height(allBlocks);
audit.trial_rows = height(allTrials);
audit.created = result.created;
write_json(fullfile(outputRoot, 'group_stats_manifest.json'), audit);
end

function cfg = load_config(configFile)
cfg = jsondecode(fileread(configFile));
required = {'raw_root_env','raw_root','processed_root_env','processed_root', ...
    'results_root_env','results_root','fieldtrip_root_env','fieldtrip_root', ...
    'megconnectome_root_env','megconnectome_root','results_subdir', ...
    'subject_manifest_relative','scan_ids','frame_interval_seconds', ...
    'min_window_coverage_fraction','bootstrap_resamples','bootstrap_seed', ...
    'confidence_level','alpha'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i}), error('MEG:WorkflowGroupConfig', 'Missing field %s.', required{i}); end
end
cfg.raw_root = resolve_root(cfg, 'raw_root_env', 'raw_root');
cfg.processed_root = resolve_root(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_root(cfg, 'results_root_env', 'results_root');
cfg.fieldtrip_root = resolve_root(cfg, 'fieldtrip_root_env', 'fieldtrip_root');
cfg.megconnectome_root = resolve_root(cfg, 'megconnectome_root_env', 'megconnectome_root');
end

function subjects = read_subject_manifest(filename)
if ~isfile(filename), error('MEG:WorkflowSubjectManifest', 'Subject manifest not found: %s', filename); end
fid = fopen(filename, 'r');
cleanup = onCleanup(@() fclose(fid));
data = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
subjects = data{1};
subjects = subjects(~cellfun(@isempty, subjects));
for i = 1:numel(subjects)
    subjects{i} = strtrim(subjects{i});
    if isempty(regexp(subjects{i}, '^[0-9]+$', 'once'))
        error('MEG:WorkflowSubjectManifest', 'Invalid subject identifier: %s', subjects{i});
    end
end
subjects = unique(subjects, 'stable');
clear cleanup
end

function validate_subject_table(tbl, subjectId, scanId, kind)
required = {'subject_id','scan_id','condition_code','valid_window', ...
    'coverage_fraction','PQ_mean','PS_mean'};
for i = 1:numel(required)
    if ~ismember(required{i}, tbl.Properties.VariableNames)
        error('MEG:WorkflowTableColumns', '%s table lacks %s.', kind, required{i});
    end
end
ids = identifiers_to_cellstr(tbl.subject_id);
scans = identifiers_to_cellstr(tbl.scan_id);
if any(~strcmp(ids, subjectId)) || any(~strcmp(scans, scanId))
    error('MEG:WorkflowTableIdentity', '%s table identity mismatch.', kind);
end
if any(~ismember(tbl.condition_code, [1 2 3]))
    error('MEG:WorkflowTableCondition', 'Unexpected condition code in %s table.', kind);
end
end

function rows = aggregate_scan_means(tbl, subjects, scanIds)
rows = repmat(empty_scan_row(), numel(subjects) * numel(scanIds) * 3, 1);
% Convert the identity columns once.  Repeating num2str/cellfun over the
% complete trial table inside the nested condition loop is prohibitively slow
% on MATLAB R2018b and a network filesystem.
ids = identifiers_to_cellstr(tbl.subject_id);
scans = identifiers_to_cellstr(tbl.scan_id);
count = 0;
for iSubject = 1:numel(subjects)
    for iScan = 1:numel(scanIds)
        for code = 1:3
            count = count + 1;
            subjectId = subjects{iSubject}; scanId = scanIds{iScan};
            ix = strcmp(ids, subjectId) & strcmp(scans, scanId) & ...
                tbl.condition_code == code & tbl.valid_window;
            rows(count).subject_id = subjectId;
            rows(count).scan_id = scanId;
            rows(count).condition_code = code;
            rows(count).condition = condition_name(code);
            rows(count).n_windows = nnz(ix);
            if rows(count).n_windows > 0
                rows(count).coverage_mean = mean(tbl.coverage_fraction(ix));
                rows(count).PQ_mean = mean(tbl.PQ_mean(ix));
                rows(count).PQ_median = median(tbl.PQ_mean(ix));
                rows(count).PS_mean = mean(tbl.PS_mean(ix));
                rows(count).PS_median = median(tbl.PS_mean(ix));
                if ismember('PQ_curl_mean', tbl.Properties.VariableNames)
                    rows(count).PQ_curl_mean = mean(tbl.PQ_curl_mean(ix));
                    rows(count).PS_curl_mean = mean(tbl.PS_curl_mean(ix));
                end
            end
        end
    end
end
rows = rows(1:count);
end

function rows = aggregate_participant_means(scanRows, subjects, scanIds)
if isstruct(scanRows)
    scanRows = struct2table(scanRows);
end
rows = repmat(empty_participant_row(), numel(subjects) * 3, 1);
ids = identifiers_to_cellstr(scanRows.subject_id);
count = 0;
for iSubject = 1:numel(subjects)
    for code = 1:3
        count = count + 1;
        rows(count).subject_id = subjects{iSubject};
        rows(count).condition_code = code;
        rows(count).condition = condition_name(code);
        ix = strcmp(ids, subjects{iSubject}) & ...
            scanRows.condition_code == code;
        selected = scanRows(ix, :);
        rows(count).n_scans_PQ = nnz(isfinite(selected.PQ_mean));
        rows(count).n_scans_PS = nnz(isfinite(selected.PS_mean));
        rows(count).n_windows = sum(selected.n_windows);
        rows(count).PQ_mean = safe_mean(selected.PQ_mean);
        rows(count).PQ_median = safe_median(selected.PQ_mean);
        rows(count).PS_mean = safe_mean(selected.PS_mean);
        rows(count).PS_median = safe_median(selected.PS_mean);
        rows(count).coverage_mean = safe_mean(selected.coverage_mean);
        rows(count).complete_PQ = rows(count).n_scans_PQ == numel(scanIds);
        rows(count).complete_PS = rows(count).n_scans_PS == numel(scanIds);
    end
end
rows = rows(1:count);
end

function stats = analyze_block_level(tbl, subjects, cfg)
stats = struct;
stats.pairwise = repmat(empty_pairwise_row(), 6, 1);
stats.descriptive = repmat(empty_descriptive_row(), 0, 1);
metrics = {'PQ','PS'};
pairCodes = [1 2; 1 3; 2 3];
pairNames = {'0-back_vs_2-back'; '0-back_vs_fixation'; '2-back_vs_fixation'};
pairCount = 0;
for iMetric = 1:2
    metricName = metrics{iMetric};
    mat = participant_matrix(tbl, subjects, metricName);
    complete = all(isfinite(mat), 2);
    x = mat(complete, :);
    for j = 1:3
        pairCount = pairCount + 1;
        a = x(:, pairCodes(j,1)); b = x(:, pairCodes(j,2));
        stats.pairwise(pairCount) = paired_test_row('block', metricName, ...
            pairNames{j}, a, b, cfg, 100*iMetric+j);
    end
    stats.descriptive(end+1) = descriptive_row('block', metricName, ...
        '0-back', x(:,1)); %#ok<AGROW>
    stats.descriptive(end+1) = descriptive_row('block', metricName, ...
        '2-back', x(:,2)); %#ok<AGROW>
    stats.descriptive(end+1) = descriptive_row('block', metricName, ...
        'fixation', x(:,3)); %#ok<AGROW>
end
stats.pairwise = apply_holm(stats.pairwise);
end

function mat = participant_matrix(tbl, subjects, metricName)
if isstruct(tbl)
    tbl = struct2table(tbl);
end
mat = NaN(numel(subjects), 3);
field = [metricName '_mean'];
ids = identifiers_to_cellstr(tbl.subject_id);
for i = 1:numel(subjects)
    for code = 1:3
        ix = strcmp(ids, subjects{i}) & ...
            tbl.condition_code == code;
        if any(ix)
            value = tbl.(field)(find(ix, 1, 'first'));
            if (strcmp(metricName, 'PQ') && tbl.complete_PQ(find(ix,1,'first'))) || ...
                    (strcmp(metricName, 'PS') && tbl.complete_PS(find(ix,1,'first')))
                mat(i, code) = value;
            end
        end
    end
end
end

function row = paired_test_row(level, metric, contrast, a, b, cfg, seedOffset)
row = empty_pairwise_row();
row.analysis_level = level;
row.metric = metric;
row.contrast = contrast;
finite = isfinite(a) & isfinite(b);
a = a(finite); b = b(finite);
d = a - b;
row.n_subjects = numel(d);
row.mean_a = safe_mean(a); row.mean_b = safe_mean(b);
row.median_a = safe_median(a); row.median_b = safe_median(b);
row.mean_difference_a_minus_b = safe_mean(d);
row.median_difference_a_minus_b = safe_median(d);
nonzero = d ~= 0;
if any(nonzero)
    row.rank_biserial = (nnz(d > 0) - nnz(d < 0)) / nnz(nonzero);
else
    row.rank_biserial = 0;
end
if numel(d) >= 1
    try
        [p, h, rs] = signrank(a, b);
        row.p_raw = p; row.reject_uncorrected = logical(h);
        if isstruct(rs) && isfield(rs, 'zval'), row.z_value = rs.zval; end
        if isstruct(rs) && isfield(rs, 'signedrank'), row.signed_rank = rs.signedrank; end
    catch
        row.p_raw = NaN;
    end
end
[row.bootstrap_ci_low, row.bootstrap_ci_high] = bootstrap_difference_ci(d, ...
    cfg.bootstrap_resamples, cfg.bootstrap_seed + seedOffset, cfg.confidence_level);
row.alpha = cfg.alpha;
end

function rows = apply_holm(rows)
p = [rows.p_raw]';
finite = isfinite(p);
adj = NaN(size(p));
if any(finite)
    idx = find(finite); [sorted, order] = sort(p(finite));
    m = numel(sorted); running = 0;
    for i = 1:m
        value = min(1, (m - i + 1) * sorted(i));
        running = max(running, value);
        adj(idx(order(i))) = running;
    end
end
for i = 1:numel(rows)
    rows(i).p_holm = adj(i);
    rows(i).reject_holm = isfinite(adj(i)) && adj(i) < rows(i).alpha;
end
end

function [lo, hi] = bootstrap_difference_ci(d, B, seed, confidence)
d = d(isfinite(d));
if numel(d) < 2 || B < 1
    lo = NaN; hi = NaN; return
end
previous = rng;
cleanup = onCleanup(@() rng(previous));
rng(seed, 'twister');
sampleIndex = randi(numel(d), numel(d), B);
boot = mean(d(sampleIndex), 1);
tail = (1 - confidence) / 2 * 100;
lo = prctile(boot, tail); hi = prctile(boot, 100 - tail);
clear cleanup
end

function row = descriptive_row(level, metric, condition, x)
row = empty_descriptive_row();
row.analysis_level = level; row.metric = metric; row.condition = condition;
x = x(isfinite(x)); row.n_subjects = numel(x);
row.mean = safe_mean(x); row.median = safe_median(x);
if numel(x) > 1, row.sd = std(x); else, row.sd = NaN; end
if ~isempty(x), row.q25 = prctile(x,25); row.q75 = prctile(x,75); end
end

function row = empty_scan_row()
row = struct('subject_id','','scan_id','','condition_code',NaN,'condition','', ...
    'n_windows',0,'coverage_mean',NaN,'PQ_mean',NaN,'PQ_median',NaN, ...
    'PS_mean',NaN,'PS_median',NaN,'PQ_curl_mean',NaN,'PS_curl_mean',NaN);
end

function row = empty_participant_row()
row = struct('subject_id','','condition_code',NaN,'condition','', ...
    'n_scans_PQ',0,'n_scans_PS',0,'n_windows',0,'coverage_mean',NaN, ...
    'PQ_mean',NaN,'PQ_median',NaN,'PS_mean',NaN,'PS_median',NaN, ...
    'complete_PQ',false,'complete_PS',false);
end

function row = empty_pairwise_row()
row = struct('analysis_level','','metric','','contrast','','n_subjects',0, ...
    'mean_a',NaN,'mean_b',NaN,'median_a',NaN,'median_b',NaN, ...
    'mean_difference_a_minus_b',NaN,'median_difference_a_minus_b',NaN, ...
    'rank_biserial',NaN,'signed_rank',NaN,'z_value',NaN,'p_raw',NaN, ...
    'p_holm',NaN,'reject_uncorrected',false,'reject_holm',false, ...
    'bootstrap_ci_low',NaN,'bootstrap_ci_high',NaN,'alpha',NaN);
end

function row = empty_descriptive_row()
row = struct('analysis_level','','metric','','condition','','n_subjects',0, ...
    'mean',NaN,'median',NaN,'sd',NaN,'q25',NaN,'q75',NaN);
end

function row = failure_row(subjectId, scanId, reason)
row = struct('subject_id',subjectId,'scan_id',scanId,'reason',reason);
end

function name = condition_name(code)
if code == 1, name = '0-back'; elseif code == 2, name = '2-back';
elseif code == 3, name = 'fixation'; else, name = 'unknown'; end
end

function ids = identifiers_to_cellstr(x)
if isnumeric(x), ids = arrayfun(@(v) num2str(v, '%.0f'), x(:), 'UniformOutput', false);
elseif iscell(x), ids = cellfun(@char, x(:), 'UniformOutput', false);
else, ids = cellstr(x); end
end

function scanIds = normalize_scan_ids(scanIds)
if isnumeric(scanIds), scanIds = arrayfun(@num2str, scanIds(:)', 'UniformOutput', false);
elseif ischar(scanIds), scanIds = {scanIds}; elseif isstring(scanIds), scanIds = cellstr(scanIds(:)'); end
end

function root = resolve_root(cfg, envField, valueField)
root = getenv(char(cfg.(envField))); if isempty(root), root = char(cfg.(valueField)); end
end

function value = safe_mean(x)
x = x(isfinite(x)); if isempty(x), value = NaN; else, value = mean(x); end
end
function value = safe_median(x)
x = x(isfinite(x)); if isempty(x), value = NaN; else, value = median(x); end
end

function write_struct_table(rows, filename)
if isempty(rows)
    return
end
rows = rows(:);
writetable(struct2table(rows), filename);
end

function write_json(filename, value)
fid = fopen(filename, 'w'); if fid < 0, error('MEG:ManifestWrite', 'Cannot open %s.', filename); end
cleanup = onCleanup(@() fclose(fid)); fwrite(fid, jsonencode(value), 'char'); clear cleanup
end

function rows = describe_conditions(tbl,subjects,level)
rows = repmat(empty_descriptive_row(),0,1);
for name={'PQ','PS'}
    x=participant_matrix(tbl,subjects,name{1});
    x=x(all(isfinite(x),2),:);
    for code=1:3
        rows(end+1)=descriptive_row(level,name{1},condition_name(code),x(:,code)); %#ok<AGROW>
    end
end
end
