function summary = meg_group_contrast_maps(configFile)
%MEG_GROUP_CONTRAST_MAPS Aggregate participant-level copes into d maps.
%
% Each participant contributes one equal-weight mean of scan 6 and scan 7.
% The group estimand follows the Project-01 level2_analysis_brainimg.m
% definition exactly: sample mean divided by sample SD plus eps.  No
% significance or FDR map is generated here.

if nargin ~= 1
    error('MEG:WorkflowGroupArguments', ...
        'Usage: meg_group_contrast_maps(configFile)');
end
cfg = load_config(configFile);
resultsRoot = resolve_root(cfg, 'results_root_env', 'results_root');
groupRoot = fullfile(resultsRoot, char(cfg.results_subdir), char(cfg.group_subdir));
ensure_directory(groupRoot);
outputFile = fullfile(groupRoot, 'cohens_d_level2.mat');
manifestFile = fullfile(groupRoot, 'group_level2_manifest.json');
if is_complete_group(outputFile, manifestFile)
    loaded = load(outputFile, 'group_result');
    if isfield(loaded, 'group_result')
        summary = loaded.group_result.summary;
        return
    end
end

subjectFile = resolve_subject_file(configFile, cfg);
if ~isfile(subjectFile)
    error('MEG:WorkflowSubjectManifest', 'Subject manifest not found: %s', subjectFile);
end
subjectIds = read_subject_ids(subjectFile);
if isfield(cfg,'expected_subjects') && numel(subjectIds) ~= double(cfg.expected_subjects)
    error('MEG:WorkflowSubjectCount', 'Expected %g subjects, found %d.', ...
        double(cfg.expected_subjects), numel(subjectIds));
end

% The registered scan estimator is an equal-weight mean of scan 6 and scan 7.
% When the user elects to exclude low-coverage scans, retain that estimand by
% including only participants with a complete pair rather than silently
% replacing a missing scan with a single-scan estimate.
includeCompleteOnly = false;
if isfield(cfg, 'group_include_complete_participants_only')
    includeCompleteOnly = logical(cfg.group_include_complete_participants_only);
end
if includeCompleteOnly
    inclusionPolicy = 'complete_scan6_scan7_pairs_only';
else
    inclusionPolicy = 'all_subjects_required';
end

copeStack = [];
V_mask = [];
maskIndex = [];
patternNames = [];
contrastNames = [];
participantFiles = cell(0,1);
includedSubjectIds = cell(0,1);
participantCopes = cell(0,1);
excludedSubjects = repmat(struct('subject_id','', 'reason','', ...
    'summary_status','', 'failed_scans',{{}}, 'error_identifiers',{{}}), 0, 1);
for iSubject = 1:numel(subjectIds)
    subjectId = subjectIds{iSubject};
    participantFile = fullfile(resultsRoot, char(cfg.results_subdir), subjectId, ...
        char(cfg.participant_subdir), 'level1_participant.mat');
    subjectSummaryFile = fullfile(resultsRoot, char(cfg.results_subdir), subjectId, ...
        'subject_summary.json');
    [subjectSummary, subjectStatus] = read_subject_summary(subjectSummaryFile);
    if includeCompleteOnly && (~strcmpi(subjectStatus, 'complete') || ...
            ~isfile(participantFile))
        if ~strcmpi(subjectStatus, 'complete')
            reason = 'subject_summary_incomplete';
        else
            reason = 'participant_level1_missing';
        end
        excludedSubjects(end+1) = make_exclusion_record(subjectId, reason, ...
            subjectStatus, subjectSummary); %#ok<AGROW>
        continue
    end
    if ~isfile(participantFile)
        error('MEG:WorkflowParticipantMissing', ...
            'Missing participant Level1 file (%d/%d): %s', ...
            iSubject, numel(subjectIds), participantFile);
    end
    loaded = load(participantFile, 'participant');
    if ~isfield(loaded, 'participant')
        error('MEG:WorkflowParticipantVariable', 'Missing participant variable: %s', participantFile);
    end
    p = loaded.participant;
    validate_participant(p, subjectId, cfg);
    if isempty(V_mask)
        V_mask = logical(p.V_mask);
        maskIndex = p.mask_index(:);
        patternNames = p.pattern_names;
        contrastNames = {p.contrasts.name};
    else
        stackSize = [numel(includedSubjectIds), nnz(V_mask), ...
            numel(contrastNames), 4];
        if ~isequal(size(p.V_mask), size(V_mask)) || ~isequal(p.V_mask, V_mask) || ...
                ~isequal(p.mask_index(:), maskIndex) || ...
                ~isequal(size(p.cope), stackSize(2:4))
            error('MEG:WorkflowParticipantShape', ...
                'Participant shape/mask mismatch for %s.', subjectId);
        end
    end
    includedSubjectIds{end+1,1} = subjectId; %#ok<AGROW>
    participantFiles{end+1,1} = participantFile; %#ok<AGROW>
    participantCopes{end+1,1} = p.cope; %#ok<AGROW>
end

if isempty(includedSubjectIds)
    error('MEG:WorkflowGroupNoSubjects', ...
        'No complete scan-6/scan-7 participant results are available.');
end
[V, C, Patt] = size(participantCopes{1});
copeStack = nan(numel(participantCopes), V, C, Patt);
for iParticipant = 1:numel(participantCopes)
    copeStack(iParticipant,:,:,:) = participantCopes{iParticipant};
end

% Compute finite-value support explicitly so an incomplete spatial vertex is
% visible in n_subjects_map rather than silently treated as zero.
[groupMeanVector, groupSdVector, nSubjectsMap] = summarize_subject_copes(copeStack);
if isfield(cfg, 'min_subjects_for_map') && ...
        any(nSubjectsMap(:) < double(cfg.min_subjects_for_map))
    error('MEG:WorkflowGroupSupport', ...
        'At least one masked vertex has fewer than %g supporting subjects.', ...
        double(cfg.min_subjects_for_map));
end
cohensDVector = groupMeanVector ./ (groupSdVector + eps);

groupMean = vector_to_map(groupMeanVector, V_mask, maskIndex);
groupSd = vector_to_map(groupSdVector, V_mask, maskIndex);
cohens_d = vector_to_map(cohensDVector, V_mask, maskIndex);
n_subjects_map = vector_to_map(nSubjectsMap, V_mask, maskIndex);

Comparison = contrastNames;
PatternNames = patternNames;
group_result = struct;
group_result.schema_version = 1;
group_result.status = 'complete';
group_result.workflow = char(cfg.workflow);
group_result.subject_manifest = subjectFile;
group_result.subject_manifest_count = numel(subjectIds);
group_result.subject_count = numel(includedSubjectIds);
group_result.subject_ids = includedSubjectIds(:)';
group_result.included_subject_ids = includedSubjectIds(:)';
group_result.all_subject_ids = subjectIds;
group_result.excluded_subjects = excludedSubjects;
group_result.inclusion_policy = inclusionPolicy;
group_result.pattern_names = PatternNames;
group_result.comparisons = Comparison;
group_result.V_mask = V_mask;
group_result.mask_index = maskIndex;
group_result.cohens_d = cohens_d;
group_result.group_mean = groupMean;
group_result.group_sd = groupSd;
group_result.n_subjects_map = n_subjects_map;
group_result.cohens_d_vector = cohensDVector;
group_result.group_mean_vector = groupMeanVector;
group_result.group_sd_vector = groupSdVector;
group_result.n_subjects_vector = nSubjectsMap;
group_result.eps = eps;
group_result.level2_formula = 'group_mean ./ (group_sd + eps), sample SD';
group_result.significance_map = false;
group_result.created = datestr(now,30);
summary = struct('workflow',char(cfg.workflow), ...
    'status','complete', 'subject_count',numel(includedSubjectIds), ...
    'subject_manifest_count',numel(subjectIds), ...
    'excluded_subject_count',numel(excludedSubjects), ...
    'excluded_scan_count',count_excluded_scans(excludedSubjects), ...
    'inclusion_policy',inclusionPolicy, ...
    'contrast_count',numel(Comparison), 'pattern_count',numel(PatternNames), ...
    'map_size',size(V_mask), 'vertex_count',nnz(V_mask), ...
    'min_subjects_map',min(nSubjectsMap(:)), ...
    'max_subjects_map',max(nSubjectsMap(:)), ...
    'zero_variance_vertices',nnz(isfinite(groupSdVector) & groupSdVector == 0), ...
    'd_formula_max_abs_error',max(abs(cohensDVector(:) - ...
        groupMeanVector(:) ./ (groupSdVector(:) + eps))), ...
    'created',group_result.created);
group_result.summary = summary;

% Keep the top-level variable names stable for downstream MATLAB/ Python
% readers, matching the names documented in the experiment contract.
group_mean = groupMean; %#ok<NASGU>
group_sd = groupSd; %#ok<NASGU>
save(outputFile, 'group_result', 'cohens_d', 'group_mean', 'group_sd', ...
    'n_subjects_map', 'Comparison', 'PatternNames', 'V_mask', 'maskIndex', ...
    'groupMeanVector', 'groupSdVector', 'cohensDVector', 'nSubjectsMap', ...
    '-v7.3');
manifest = struct('schema_version',1, 'status','complete', ...
    'workflow',char(cfg.workflow), 'subject_manifest',subjectFile, ...
    'subject_manifest_count',numel(subjectIds), ...
    'subject_count',numel(includedSubjectIds), ...
    'subject_ids',{includedSubjectIds(:)'}, ...
    'included_subject_ids',{includedSubjectIds(:)'}, ...
    'all_subject_ids',{subjectIds}, 'excluded_subjects',excludedSubjects, ...
    'inclusion_policy',inclusionPolicy, ...
    'participant_files',{participantFiles}, 'comparisons',{Comparison}, ...
    'pattern_names',{PatternNames}, 'mask_size',size(V_mask), ...
    'vertex_count',nnz(V_mask), 'cohens_d_dimensions',size(cohens_d), ...
    'group_mean_dimensions',size(groupMean), 'group_sd_dimensions',size(groupSd), ...
    'n_subjects_map_dimensions',size(n_subjects_map), ...
    'min_subjects_map',min(nSubjectsMap(:)), 'max_subjects_map',max(nSubjectsMap(:)), ...
    'zero_variance_vertices',nnz(isfinite(groupSdVector) & groupSdVector == 0), ...
    'level2_formula',group_result.level2_formula, 'significance_map',false, ...
    'created',group_result.created);
write_json(manifestFile, manifest);
end

function cfg = load_config(configFile)
if ~isfile(configFile), error('MEG:WorkflowGroupConfig', 'Missing config: %s', configFile); end
cfg = jsondecode(fileread(configFile));
required = {'workflow','results_root_env','results_root','results_subdir', ...
    'participant_subdir','group_subdir','subject_manifest_relative'};
for i = 1:numel(required)
    if ~isfield(cfg,required{i}), error('MEG:WorkflowGroupConfigField', ...
            'Missing config field %s.', required{i}); end
end
end

function root = resolve_root(cfg, envField, fallbackField)
root = getenv(char(cfg.(envField)));
if isempty(root), root = char(cfg.(fallbackField)); end
end

function subjectFile = resolve_subject_file(configFile, cfg)
override = getenv('COF_MEG_SUBJECT_MANIFEST');
if ~isempty(override), subjectFile = override; return; end
projectRoot = fileparts(fileparts(char(configFile)));
subjectFile = meg_subject_manifest(cfg, projectRoot);
end

function ids = read_subject_ids(filename)
lines = strsplit(strtrim(fileread(filename)), {'\r','\n'});
lines = lines(~cellfun(@isempty, lines));
ids = cell(size(lines));
for i = 1:numel(lines)
    ids{i} = char(strtrim(lines{i}));
    if isempty(regexp(ids{i}, '^[A-Za-z0-9_-]+$', 'once'))
        error('MEG:WorkflowSubjectIdentifier', 'Invalid subject entry: %s', ids{i});
    end
end
ids = ids(:)';
if numel(unique(ids)) ~= numel(ids)
    error('MEG:WorkflowSubjectIdentifier', 'Subject manifest contains duplicate identifiers.');
end
end

function validate_participant(p, subjectId, cfg)
required = {'status','subject_id','scan_ids','scan_aggregation','cope', ...
    'V_mask','mask_index','pattern_names','contrasts'};
for i = 1:numel(required)
    if ~isfield(p,required{i}), error('MEG:WorkflowParticipantField', ...
            'Participant %s lacks %s.', subjectId, required{i}); end
end
if ~strcmpi(char(p.status),'complete') || ~strcmp(char(p.subject_id),subjectId)
    error('MEG:WorkflowParticipantStatus', 'Participant status/id mismatch for %s.', subjectId);
end
if ~strcmp(char(p.scan_aggregation),'equal_weight_mean_scan6_scan7') || ...
        numel(p.scan_ids) ~= 2 || ~all(ismember({'6','7'}, normalize_scan_ids(p.scan_ids)))
    error('MEG:WorkflowParticipantScans', 'Participant %s is not a scan-6/7 equal-weight result.', subjectId);
end
if ~isequal(size(p.cope), [nnz(p.V_mask), numel(cfg.contrasts), 4])
    error('MEG:WorkflowParticipantShape', 'Unexpected cope shape for %s.', subjectId);
end
if any(~isfinite(p.cope(:)))
    error('MEG:WorkflowParticipantFinite', 'Nonfinite participant cope for %s.', subjectId);
end
end

function [groupMean, groupSd, nMap] = summarize_subject_copes(copeStack)
[N,V,C,Patt] = size(copeStack);
groupMean = nan(V,C,Patt);
groupSd = nan(V,C,Patt);
nMap = zeros(V,C,Patt);
for c = 1:C
    for p = 1:Patt
        for v = 1:V
            values = copeStack(:,v,c,p);
            keep = isfinite(values);
            n = nnz(keep);
            nMap(v,c,p) = n;
            if n > 0
                values = values(keep);
                groupMean(v,c,p) = mean(values);
                if n > 1, groupSd(v,c,p) = std(values,0); end
            end
        end
    end
end
end

function mapped = vector_to_map(vector, V_mask, maskIndex)
mapSize = size(V_mask);
[~,C,Patt] = size(vector);
mapped = nan(mapSize(1), mapSize(2), C, Patt);
for c = 1:C
    for p = 1:Patt
        tmp = nan(numel(V_mask),1);
        tmp(maskIndex) = vector(:,c,p);
        mapped(:,:,c,p) = reshape(tmp, mapSize);
    end
end
end

function ids = normalize_scan_ids(ids)
if ischar(ids), ids = {ids}; elseif isstring(ids), ids = cellstr(ids(:)'); end
for i = 1:numel(ids), ids{i} = char(string(ids{i})); end
end

function [summary, status] = read_subject_summary(filename)
summary = struct;
status = 'missing';
if ~isfile(filename), return; end
try
    summary = jsondecode(fileread(filename));
    if isfield(summary, 'status')
        status = char(summary.status);
    else
        status = 'invalid';
    end
catch
    status = 'invalid';
end
end

function record = make_exclusion_record(subjectId, reason, summaryStatus, summary)
record = struct('subject_id',subjectId, 'reason',reason, ...
    'summary_status',summaryStatus, 'failed_scans',{{}}, ...
    'error_identifiers',{{}});
if ~isstruct(summary) || ~isfield(summary, 'scans'), return; end
failedScans = cell(0,1);
errorIdentifiers = cell(0,1);
for iScan = 1:numel(summary.scans)
    scan = summary.scans(iScan);
    if isfield(scan,'status') && ~strcmpi(char(scan.status),'complete')
        if isfield(scan,'scan_id')
            failedScans{end+1,1} = char(scan.scan_id); %#ok<AGROW>
        end
        if isfield(scan,'error_identifier') && ~isempty(scan.error_identifier)
            errorIdentifiers{end+1,1} = char(scan.error_identifier); %#ok<AGROW>
        end
    end
end
record.failed_scans = failedScans(:)';
record.error_identifiers = errorIdentifiers(:)';
end

function n = count_excluded_scans(records)
n = 0;
for i = 1:numel(records)
    n = n + numel(records(i).failed_scans);
end
end

function tf = is_complete_group(outputFile, manifestFile)
tf = false;
if ~isfile(outputFile) || ~isfile(manifestFile), return; end
manifest = jsondecode(fileread(manifestFile));
tf = isfield(manifest,'status') && strcmpi(char(manifest.status),'complete');
end

function ensure_directory(directory)
if ~isfolder(directory), mkdir(directory); end
end

function write_json(filename, value)
fid = fopen(filename,'w');
if fid < 0, error('MEG:WorkflowGroupManifest', 'Cannot open %s.', filename); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value), 'char');
clear cleanup
end
