function summary = meg_native_contrast_maps(configFile, subjectId, scanIds)
%MEG_NATIVE_CONTRAST_MAPS Run meg_native_contrast_maps level-1 regression for one subject.
%
% The parent P04-EXP-003 LCMV/COF derivatives are treated as immutable input.
% This entry point materialises D/C maps with the registered P04 parameters,
% decodes official task/fixation blocks, performs a segment-aware GLS-AR(1)
% regression for each scan, and averages scan 6/7 copes with equal weight.
% Subject-level arrays stay under the configured (normally HPC) results root.

if nargin < 2 || nargin > 3
    error('MEG:WorkflowArguments', ...
        'Usage: meg_native_contrast_maps(configFile, subjectId, [scanIds])');
end
cfg = load_config(configFile);
codeRoot = fileparts(mfilename('fullpath'));
addpath(codeRoot);
addpath(cfg.fieldtrip_root);
ft_defaults;
addpath(fullfile(cfg.megconnectome_root, 'trial_functions'));

subjectId = validate_identifier(subjectId, 'subjectId');
if nargin < 3 || isempty(scanIds), scanIds = cfg.scan_ids; end
scanIds = normalize_scan_ids(scanIds);
if numel(scanIds) ~= 2 || ~all(ismember({'6','7'}, scanIds))
    error('MEG:WorkflowScans', ...
        'meg_native_contrast_maps requires exactly scan 6 and scan 7 for equal weighting.');
end

resultsRoot = resolve_root(cfg, 'results_root_env', 'results_root');
subjectRoot = fullfile(resultsRoot, char(cfg.results_subdir), subjectId);
ensure_directory(subjectRoot);
metricOptions = default_metric_options(cfg);
glsOptions = struct('include_unlabeled_intervals', ...
    logical(cfg.include_unlabeled_intervals_as_intercept), 'min_observations', 3);

scanSummary = repmat(struct('scan_id','', 'status','', ...
    'source_directory','', 'cof_directory','', 'maps_directory','', ...
    'level1_file','', 'n_blocks',0, 'n_map_intervals',0, ...
    'frame_count',0, 'condition_counts',[], 'error_identifier','', ...
    'error_message',''), 1, numel(scanIds));
scanResults = cell(1, numel(scanIds));

for iScan = 1:numel(scanIds)
    scanId = scanIds{iScan};
    scanSummary(iScan).scan_id = scanId;
    scanRoot = fullfile(subjectRoot, ['scan-' scanId]);
    ensure_directory(scanRoot);
    [sourceDir, cofDir] = source_directories(cfg, subjectId, scanId);
    mapsDir = fullfile(scanRoot, char(cfg.map_metrics_subdir));
    level1Dir = fullfile(scanRoot, char(cfg.level1_subdir));
    level1File = fullfile(level1Dir, 'level1_scan.mat');
    scanSummary(iScan).source_directory = sourceDir;
    scanSummary(iScan).cof_directory = cofDir;
    scanSummary(iScan).maps_directory = mapsDir;
    scanSummary(iScan).level1_file = level1File;
    try
        require_complete_cof(cofDir, sourceDir, subjectId, scanId);
        if ~is_complete_maps(mapsDir)
            move_incomplete_directory_aside(mapsDir);
            meg_native_motif_metrics(cofDir, mapsDir, [1 Inf], metricOptions);
        end

        [V_mask, maskInfo] = load_cof_mask(cofDir);
        rawFile = resolve_raw_file(cfg, subjectId, scanId);
        metadataFile = fullfile(sourceDir, 'session_metadata.mat');
        if ~isfile(metadataFile)
            error('MEG:WorkflowMetadata', 'Missing session metadata: %s', metadataFile);
        end
        loadedMetadata = load(metadataFile, 'metadata');
        if ~isfield(loadedMetadata, 'metadata') || ...
                ~isfield(loadedMetadata.metadata, 'session_input_info')
            error('MEG:WorkflowMetadata', ...
                'session_metadata.mat lacks session_input_info: %s', metadataFile);
        end
        sessionInfo = loadedMetadata.metadata.session_input_info;
        decoded = meg_decode_events(rawFile, cfg.megconnectome_root, 0, 0);
        [blockTimes, blockQC] = make_block_times(decoded, sessionInfo, cfg);
        blockQC = add_map_coverage_qc(blockQC, mapsDir, blockTimes, ...
            cfg.frame_interval_seconds, cfg.min_block_coverage_fraction);

        if ~is_complete_level1(level1File)
            ensure_directory(level1Dir);
            result = meg_fit_streaming_gls(mapsDir, V_mask, ...
                blockTimes, cfg.frame_interval_seconds, cfg.contrasts, glsOptions);
            result.workflow = char(cfg.workflow);
            result.subject_id = subjectId;
            result.scan_id = scanId;
            result.parent_derivative_experiment = char(cfg.parent_workflow);
            result.source_workflow = char(cfg.source_workflow);
            result.source_cof_directory = cofDir;
            result.maps_directory = mapsDir;
            result.block_times = blockTimes;
            result.block_qc = blockQC;
            result.mask_info = maskInfo;
            save(level1File, 'result', '-v7.3');
            level1Manifest = make_level1_manifest(result, cfg, subjectId, scanId, ...
                sourceDir, cofDir, mapsDir, metadataFile, rawFile, blockQC);
            write_json(fullfile(level1Dir, 'level1_scan_manifest.json'), level1Manifest);
        else
            loadedResult = load(level1File, 'result');
            if ~isfield(loadedResult, 'result')
                error('MEG:WorkflowLevel1File', 'Missing result variable: %s', level1File);
            end
            result = loadedResult.result;
            validate_level1_result(result, V_mask, cfg, subjectId, scanId);
            % Recompute coverage/segment diagnostics on every resume so that
            % an older complete Level-1 file is upgraded to the current
            % manifest contract without rerunning GLS.
            result.block_qc = blockQC;
            level1Manifest = make_level1_manifest(result, cfg, subjectId, scanId, ...
                sourceDir, cofDir, mapsDir, metadataFile, rawFile, blockQC);
            save(level1File, 'result', '-v7.3');
            write_json(fullfile(level1Dir, 'level1_scan_manifest.json'), level1Manifest);
        end

        scanResults{iScan} = result;
        scanSummary(iScan).status = 'complete';
        scanSummary(iScan).n_blocks = blockQC.n_blocks;
        scanSummary(iScan).n_map_intervals = sum(result.frame_count);
        scanSummary(iScan).frame_count = result.frame_count;
        scanSummary(iScan).condition_counts = result.condition_counts;
    catch ME
        scanSummary(iScan).status = 'failed';
        scanSummary(iScan).error_identifier = ME.identifier;
        scanSummary(iScan).error_message = getReport(ME, 'extended', 'hyperlinks', 'off');
        failure = struct('schema_version',1, 'status','failed', ...
            'workflow',char(cfg.workflow), 'subject_id',subjectId, ...
            'scan_id',scanId, 'error_identifier',ME.identifier, ...
            'error_message',scanSummary(iScan).error_message, ...
            'created',datestr(now,30));
        write_json(fullfile(scanRoot, 'level1_scan_manifest.json'), failure);
    end
end

complete = strcmp({scanSummary.status}, 'complete');
participantFile = fullfile(subjectRoot, char(cfg.participant_subdir), ...
    'level1_participant.mat');
if all(complete)
    if ~is_complete_participant(participantFile)
        participantDir = fileparts(participantFile);
        ensure_directory(participantDir);
        copeByScan = cat(4, scanResults{1}.cope, scanResults{2}.cope);
        cope = mean(copeByScan, 4);
        participant = struct;
        participant.schema_version = 1;
        participant.status = 'complete';
        participant.workflow = char(cfg.workflow);
        participant.subject_id = subjectId;
        participant.scan_ids = scanIds;
        participant.scan_aggregation = 'equal_weight_mean_scan6_scan7';
        participant.pattern_names = scanResults{1}.pattern_names;
        participant.contrasts = cfg.contrasts;
        participant.contrast_vector = scanResults{1}.contrast_vector;
        participant.V_mask = scanResults{1}.V_mask;
        participant.mask_index = scanResults{1}.mask_index;
        participant.cope_by_scan = copeByScan;
        participant.cope = cope;
        participant.created = datestr(now,30);
        save(participantFile, 'participant', '-v7.3');
        participantManifest = struct('schema_version',1, 'status','complete', ...
            'workflow',char(cfg.workflow), 'subject_id',subjectId, ...
            'scan_ids',{scanIds}, 'scan_aggregation',participant.scan_aggregation, ...
            'level1_scan_files',{{scanSummary.level1_file}}, ...
            'cope_dimensions',size(cope), 'created',participant.created);
        write_json(fullfile(fileparts(participantFile), 'level1_participant_manifest.json'), ...
            participantManifest);
    else
        loadedParticipant = load(participantFile, 'participant');
        if ~isfield(loadedParticipant, 'participant')
            error('MEG:WorkflowParticipantFile', ...
                'Missing participant variable: %s', participantFile);
        end
    end
end

summary = struct;
summary.schema_version = 1;
summary.status = 'complete';
summary.workflow = char(cfg.workflow);
summary.subject_id = subjectId;
summary.scan_ids = scanIds;
summary.scans = scanSummary;
summary.participant_file = participantFile;
summary.config_file = configFile;
summary.results_root = resultsRoot;
summary.created = datestr(now,30);
if ~all(complete), summary.status = 'incomplete'; end
write_json(fullfile(subjectRoot, 'subject_summary.json'), summary);
save(fullfile(subjectRoot, 'subject_summary.mat'), 'summary');
if ~strcmp(summary.status, 'complete')
    error('MEG:WorkflowSubjectFailed', ...
        'meg_native_contrast_maps failed for %s (%d of %d scans).', ...
        subjectId, nnz(~complete), numel(complete));
end
end

function cfg = load_config(configFile)
if nargin ~= 1 || ~isfile(configFile)
    error('MEG:WorkflowConfig', 'Configuration file not found: %s', char(configFile));
end
cfg = jsondecode(fileread(configFile));
required = {'workflow','parent_workflow','processed_root_env', ...
    'processed_root','results_root_env','results_root','raw_root_env','raw_root', ...
    'fieldtrip_root_env','fieldtrip_root','megconnectome_root_env', ...
    'megconnectome_root','raw_datafile_template','scan_ids', ...
    'metric_processed_subdir','source_session_scope','source_workflow', ...
    'source_cof_subdir','results_subdir', ...
    'map_metrics_subdir','level1_subdir','participant_subdir','frame_interval_seconds', ...
    'expected_blocks_per_condition_per_scan','min_block_coverage_fraction', ...
    'include_unlabeled_intervals_as_intercept','divcurl','contrasts'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i})
        error('MEG:WorkflowConfigField', 'Missing field %s.', required{i});
    end
end
cfg.processed_root = resolve_root(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_root(cfg, 'results_root_env', 'results_root');
cfg.raw_root = resolve_root(cfg, 'raw_root_env', 'raw_root');
cfg.fieldtrip_root = resolve_root(cfg, 'fieldtrip_root_env', 'fieldtrip_root');
cfg.megconnectome_root = resolve_root(cfg, 'megconnectome_root_env', 'megconnectome_root');
if cfg.frame_interval_seconds <= 0 || cfg.min_block_coverage_fraction <= 0 || ...
        cfg.min_block_coverage_fraction > 1
    error('MEG:WorkflowConfigValue', 'Invalid frame or block-coverage setting.');
end
if ~isfolder(cfg.fieldtrip_root) || ~isfolder(cfg.megconnectome_root)
    error('MEG:WorkflowDependencyRoot', ...
        'FieldTrip or megconnectome root is not available.');
end
end

function root = resolve_root(cfg, envField, fallbackField)
root = getenv(char(cfg.(envField)));
if isempty(root), root = char(cfg.(fallbackField)); end
end

function [sourceDir, cofDir] = source_directories(cfg, subjectId, scanId)
scope = lower(char(cfg.source_session_scope));
outputName = sprintf('%s_MEG_%s-%s_%s_sessionlcmv_fixedori_power_dt%gms', ...
    subjectId, scanId, char(cfg.source_workflow), scope, ...
    1000 * cfg.frame_interval_seconds);
sourceDir = fullfile(cfg.processed_root, subjectId, ...
    char(cfg.source_workflow), char(cfg.metric_processed_subdir), outputName);
cofDir = fullfile(sourceDir, char(cfg.source_cof_subdir));
end

function rawFile = resolve_raw_file(cfg, subjectId, scanId)
rawFile = char(cfg.raw_datafile_template);
rawFile = strrep(rawFile, '{raw_root}', char(cfg.raw_root));
rawFile = strrep(rawFile, '{subject_id}', subjectId);
rawFile = strrep(rawFile, '{scan_id}', scanId);
rawFile = strrep(rawFile, '{workflow}', char(cfg.source_workflow));
if ~isfile(rawFile)
    error('MEG:WorkflowRawMissing', 'Raw MEG file not found: %s', rawFile);
end
end

function require_complete_cof(cofDir, sourceDir, subjectId, scanId)
if ~isfolder(sourceDir) || ~isfolder(cofDir)
    error('MEG:WorkflowCOFDirectory', ...
        'Missing parent source/COF directory for %s scan %s.', subjectId, scanId);
end
sourceManifestFile = fullfile(sourceDir, 'session_manifest.json');
manifestFile = fullfile(cofDir, 'cof_manifest.json');
prepFile = fullfile(cofDir, 'cof_preprocessing.mat');
if ~isfile(sourceManifestFile) || ~isfile(manifestFile) || ~isfile(prepFile)
    error('MEG:WorkflowCOFFiles', 'Incomplete parent COF derivative: %s', cofDir);
end
sourceManifest = jsondecode(fileread(sourceManifestFile));
if ~isfield(sourceManifest,'status') || ~strcmpi(char(sourceManifest.status),'complete')
    error('MEG:WorkflowSourceStatus', 'Parent source manifest is not complete: %s', sourceManifestFile);
end
manifest = jsondecode(fileread(manifestFile));
if ~isfield(manifest,'status') || ~strcmpi(char(manifest.status),'complete') || ...
        ~isfield(manifest,'velocity_chunks') || isempty(manifest.velocity_chunks)
    error('MEG:WorkflowCOFStatus', 'Parent COF manifest is not complete: %s', manifestFile);
end
end

function [V_mask, maskInfo] = load_cof_mask(cofDir)
loaded = load(fullfile(cofDir, 'cof_preprocessing.mat'), 'preprocessing');
if ~isfield(loaded, 'preprocessing') || ~isfield(loaded.preprocessing, 'mask')
    error('MEG:WorkflowMask', 'COF preprocessing lacks mask: %s', cofDir);
end
mask = double(loaded.preprocessing.mask);
V_mask = conv2(mask, ones(2,2) ./ 4, 'valid');
V_mask(V_mask ~= 1) = 0;
V_mask = logical(V_mask);
if ndims(V_mask) ~= 2 || ~any(V_mask(:))
    error('MEG:WorkflowMask', 'The derived vertex mask is empty.');
end
maskInfo = struct('input_size',size(mask), 'vertex_size',size(V_mask), ...
    'vertex_count',nnz(V_mask), ...
    'construction','conv2(mask,ones(2,2)/4,''valid''); values~=1 -> 0');
end

function [blockTimes, qc] = make_block_times(decoded, sessionInfo, cfg)
if ~isfield(sessionInfo,'session_start_sample_raw')
    error('MEG:WorkflowMetadata', 'session_start_sample_raw is missing.');
end
fs = double(decoded.recording_fsample_hz);
sessionStart = double(sessionInfo.session_start_sample_raw);
blocks = decoded.blocks;
blockTimes = repmat(struct('block_index',0, 'condition_code',0, ...
    'start_sample',0, 'end_sample',0, 'start_seconds',0, 'end_seconds',0), ...
    1, numel(blocks));
for i = 1:numel(blocks)
    startSample = double(blocks(i).start_sample);
    endSample = double(blocks(i).end_sample);
    blockTimes(i).block_index = i;
    blockTimes(i).condition_code = double(blocks(i).condition_code);
    blockTimes(i).start_sample = startSample;
    blockTimes(i).end_sample = endSample;
    blockTimes(i).start_seconds = (startSample - sessionStart) ./ fs;
    blockTimes(i).end_seconds = (endSample - sessionStart + 1) ./ fs;
end
counts = zeros(1,3);
for code = 1:3, counts(code) = nnz([blockTimes.condition_code] == code); end
expected = double(cfg.expected_blocks_per_condition_per_scan);
if any(counts ~= expected)
    error('MEG:WorkflowBlockCount', ...
        'Expected %g blocks for each condition, observed [%s].', expected, num2str(counts));
end
if any([blockTimes.end_seconds] <= [blockTimes.start_seconds]) || ...
        any(~isfinite([blockTimes.start_seconds blockTimes.end_seconds]))
    error('MEG:WorkflowBlockBounds', 'Invalid decoded block bounds.');
end
qc = struct('n_blocks',numel(blockTimes), 'condition_block_counts',counts, ...
    'expected_blocks_per_condition',expected, 'recording_fsample_hz',fs, ...
    'session_start_sample_raw',sessionStart, 'min_block_coverage_fraction', ...
    cfg.min_block_coverage_fraction, 'decoder_warning_count',decoded.warning_count);
end

function qc = add_map_coverage_qc(qc, mapsDir, blockTimes, dt, minimumFraction)
manifestFile = fullfile(mapsDir, 'divcurl_manifest.json');
manifest = jsondecode(fileread(manifestFile));
allTimes = zeros(0,1);
allFrames = zeros(0,1);
for iChunk = 1:numel(manifest.chunks)
    filename = fullfile(mapsDir, char(manifest.chunks(iChunk).filename));
    loaded = load(filename, 'interval_frame_index', 'interval_time_seconds');
    allFrames = [allFrames; double(loaded.interval_frame_index(:))]; %#ok<AGROW>
    allTimes = [allTimes; double(loaded.interval_time_seconds(:))]; %#ok<AGROW>
end
[allFrames, order] = sort(allFrames);
allTimes = allTimes(order);
if isempty(allTimes)
    error('MEG:WorkflowMapCoverage', 'D/C map output contains no intervals.');
end
if any(diff(allFrames) <= 0)
    error('MEG:WorkflowFrameOrder', 'D/C interval frame indices are not strictly increasing.');
end
segmentStart = [true; diff(allFrames) ~= 1];
segmentEnd = [find(segmentStart(2:end))-1; numel(allFrames)];
segmentStartIndex = find(segmentStart);
qc.n_valid_segments = numel(segmentStartIndex);
qc.segment_frame_starts = allFrames(segmentStartIndex)';
qc.segment_frame_ends = allFrames(segmentEnd)';
qc.segment_interval_counts = (segmentEnd - segmentStartIndex + 1)';
coverageCounts = zeros(1,numel(blockTimes));
coverageExpected = zeros(1,numel(blockTimes));
coverageFraction = zeros(1,numel(blockTimes));
tol = max(1e-9, dt * 1e-4);
for iBlock = 1:numel(blockTimes)
    in = allTimes >= blockTimes(iBlock).start_seconds - tol & ...
        allTimes + dt <= blockTimes(iBlock).end_seconds + tol;
    coverageCounts(iBlock) = nnz(in);
    coverageExpected(iBlock) = max(1, round((blockTimes(iBlock).end_seconds - ...
        blockTimes(iBlock).start_seconds) ./ dt));
    coverageFraction(iBlock) = coverageCounts(iBlock) ./ coverageExpected(iBlock);
end
if any(coverageFraction < minimumFraction)
    error('MEG:WorkflowBlockCoverage', ...
        'At least one official block has map coverage below %.3f.', minimumFraction);
end
qc.map_intervals_total = numel(allTimes);
qc.map_interval_frame_start = allFrames(1);
qc.map_interval_frame_end = allFrames(end);
qc.block_interval_counts = coverageCounts;
qc.block_interval_expected = coverageExpected;
qc.block_coverage_fraction = coverageFraction;
qc.min_observed_block_coverage_fraction = min(coverageFraction);
end

function options = default_metric_options(cfg)
options = struct('filtSigma',3, 'doZscore',true, 'doBinarize',true, ...
    'nStd',2, 'sizeThreshold',10, 'saveMaps',true);
if isfield(cfg.divcurl,'filt_sigma'), options.filtSigma = cfg.divcurl.filt_sigma; end
if isfield(cfg.divcurl,'do_zscore'), options.doZscore = cfg.divcurl.do_zscore; end
if isfield(cfg.divcurl,'do_binarize'), options.doBinarize = cfg.divcurl.do_binarize; end
if isfield(cfg.divcurl,'n_std'), options.nStd = cfg.divcurl.n_std; end
if isfield(cfg.divcurl,'size_threshold'), options.sizeThreshold = cfg.divcurl.size_threshold; end
options.saveMaps = true;
end

function manifest = make_level1_manifest(result, cfg, subjectId, scanId, ...
        sourceDir, cofDir, mapsDir, metadataFile, rawFile, blockQC)
manifest = struct;
manifest.schema_version = 1;
manifest.status = 'complete';
manifest.workflow = char(cfg.workflow);
manifest.subject_id = subjectId;
manifest.scan_id = scanId;
manifest.solver = result.solver;
manifest.segment_aware = result.segment_aware;
manifest.source_directory = sourceDir;
manifest.cof_directory = cofDir;
manifest.maps_directory = mapsDir;
manifest.session_metadata_file = metadataFile;
manifest.raw_file = rawFile;
manifest.frame_interval_seconds = cfg.frame_interval_seconds;
manifest.pattern_names = result.pattern_names;
manifest.contrasts = result.contrasts;
manifest.cope_dimensions = size(result.cope);
manifest.varcope_dimensions = size(result.varcope);
manifest.tstat_dimensions = size(result.tstat);
manifest.mask_size = size(result.V_mask);
manifest.vertex_count = numel(result.mask_index);
manifest.frame_count = result.frame_count;
manifest.whitened_frame_count = result.whitened_frame_count;
manifest.condition_counts = result.condition_counts;
manifest.rho_summary = struct('min',min(result.rho_hat(:)), ...
    'max',max(result.rho_hat(:)), 'mean',mean(result.rho_hat(:)));
manifest.block_qc = blockQC;
manifest.created = datestr(now,30);
end

function validate_level1_result(result, V_mask, cfg, subjectId, scanId)
required = {'cope','varcope','tstat','beta_gls','rho_hat','V_mask', ...
    'mask_index','pattern_names','contrasts','condition_counts'};
for i = 1:numel(required)
    if ~isfield(result, required{i})
        error('MEG:WorkflowLevel1Fields', 'Level1 result lacks %s.', required{i});
    end
end
if ~isequal(size(result.V_mask), size(V_mask)) || ~isequal(result.V_mask, V_mask)
    error('MEG:WorkflowLevel1Mask', 'Level1 mask mismatch for %s scan %s.', subjectId, scanId);
end
if ~isequal(size(result.cope), [nnz(V_mask), numel(cfg.contrasts), 4])
    error('MEG:WorkflowLevel1Shape', 'Unexpected cope dimensions for %s scan %s.', subjectId, scanId);
end
end

function tf = is_complete_maps(directory)
tf = false;
manifestFile = fullfile(directory, 'divcurl_manifest.json');
progressFile = fullfile(directory, 'divcurl_progress.json');
if ~isfile(manifestFile) || ~isfile(progressFile), return; end
manifest = jsondecode(fileread(manifestFile));
progress = jsondecode(fileread(progressFile));
tf = isfield(manifest,'status') && isfield(progress,'status') && ...
    strcmpi(char(manifest.status),'complete') && strcmpi(char(progress.status),'complete') && ...
    isfield(manifest,'chunks') && ~isempty(manifest.chunks);
if tf
    for iChunk = 1:numel(manifest.chunks)
        if ~isfile(fullfile(directory, char(manifest.chunks(iChunk).filename)))
            tf = false;
            return
        end
    end
end
end

function tf = is_complete_level1(filename)
tf = false;
manifestFile = fullfile(fileparts(filename), 'level1_scan_manifest.json');
if ~isfile(filename) || ~isfile(manifestFile), return; end
manifest = jsondecode(fileread(manifestFile));
tf = isfield(manifest,'status') && strcmpi(char(manifest.status),'complete');
end

function tf = is_complete_participant(filename)
tf = false;
manifestFile = fullfile(fileparts(filename), 'level1_participant_manifest.json');
if ~isfile(filename) || ~isfile(manifestFile), return; end
manifest = jsondecode(fileread(manifestFile));
tf = isfield(manifest,'status') && strcmpi(char(manifest.status),'complete');
end

function move_incomplete_directory_aside(directory)
if ~isfolder(directory), return; end
if is_complete_maps(directory), return; end
partial = [directory '.partial'];
if isfolder(partial)
    % Preserve an earlier failed attempt while allowing an automatic retry.
    % The timestamped suffix avoids deleting or overwriting any diagnostic
    % files from a previous Slurm task.
    partial = [directory '.partial_' datestr(now,30)];
    serial = 1;
    while isfolder(partial)
        partial = sprintf('%s.partial_%s_%03d', directory, datestr(now,30), serial);
        serial = serial + 1;
    end
end
[ok,msg] = movefile(directory, partial);
if ~ok, error('MEG:WorkflowPartialMove', '%s', msg); end
end

function ensure_directory(directory)
if ~isfolder(directory), mkdir(directory); end
end

function scanIds = normalize_scan_ids(scanIds)
if isnumeric(scanIds)
    scanIds = arrayfun(@(x) num2str(x), scanIds(:)', 'UniformOutput', false);
elseif ischar(scanIds)
    scanIds = {scanIds};
elseif isstring(scanIds)
    scanIds = cellstr(scanIds(:)');
elseif ~iscell(scanIds)
    error('MEG:WorkflowScanType', 'scanIds must be numeric, char, string, or cell.');
end
for i = 1:numel(scanIds)
    scanIds{i} = validate_identifier(scanIds{i}, 'scanId');
end
end

function identifier = validate_identifier(identifier, name)
identifier = char(string(identifier));
if isempty(regexp(identifier, '^[A-Za-z0-9_-]+$', 'once'))
    error('MEG:WorkflowIdentifier', '%s contains unsupported characters.', name);
end
end

function write_json(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('MEG:WorkflowManifestWrite', 'Cannot open %s.', filename); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid, jsonencode(value), 'char');
clear cleanup
end
