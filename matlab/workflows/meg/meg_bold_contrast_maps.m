function summary = meg_bold_contrast_maps(configFile, subjectId, scanIds)
%MEG_BOLD_CONTRAST_MAPS Regress meg_bold motif maps for one participant.
%
% This function is intentionally read-only with respect to meg_bold. It
% consumes the locked exact_project1 D/C chunks and official block tables,
% runs segmented GLS-AR(1), then averages scan-6/scan-7 copes equally.

if nargin < 2 || nargin > 3
    error('MEG:WorkflowArguments', ...
        'Usage: meg_bold_contrast_maps(configFile,subjectId,[scanIds])');
end
cfg = load_config(configFile);
subjectId = validate_identifier(subjectId, 'subjectId');
if nargin < 3 || isempty(scanIds), scanIds = cfg.scan_ids; end
scanIds = normalize_scan_ids(scanIds);
if numel(scanIds) ~= 2 || ~all(ismember({'6','7'}, scanIds))
    error('MEG:WorkflowScans', 'Exactly scans 6 and 7 are required.');
end

resultsRoot = resolve_root(cfg, 'results_root_env', 'results_root');
processedRoot = resolve_root(cfg, 'processed_root_env', 'processed_root');
subjectRoot = fullfile(resultsRoot, char(cfg.results_subdir), subjectId);
ensure_directory(subjectRoot);
scanResults = cell(1,2);
scanSummary = repmat(struct('scan_id','', 'status','', 'metric_directory','', ...
    'block_table','', 'level1_file','', 'frame_count',[], ...
    'condition_counts',[], 'error_identifier','', 'error_message',''), 1, 2);

for iScan = 1:2
    scanId = scanIds{iScan};
    scanSummary(iScan).scan_id = scanId;
    metricDir = fullfile(processedRoot, subjectId, 'Wrkmem', ...
        char(cfg.derived_subdir), ['scan-' scanId], char(cfg.cof_directory), ...
        char(cfg.metric_subdir));
    blockFile = fullfile(resultsRoot, char(cfg.parent_results_subdir), ...
        subjectId, ['scan-' scanId], char(cfg.block_table_name));
    level1Dir = fullfile(subjectRoot, ['scan-' scanId], char(cfg.level1_subdir));
    level1File = fullfile(level1Dir, 'level1_scan.mat');
    scanSummary(iScan).metric_directory = metricDir;
    scanSummary(iScan).block_table = blockFile;
    scanSummary(iScan).level1_file = level1File;
    try
        validate_parent(cfg, resultsRoot, subjectId, scanId, metricDir, blockFile);
        [V_mask, metricManifest] = load_metric_mask(metricDir);
        blockTimes = read_block_times(blockFile, subjectId, scanId, cfg);
        if is_complete_level1(level1File, cfg, subjectId, scanId)
            loaded = load(level1File, 'result');
            result = loaded.result;
            validate_level1(result, V_mask, cfg, subjectId, scanId);
        else
            ensure_directory(level1Dir);
            options = struct('include_unlabeled_intervals', ...
                logical(cfg.include_unlabeled_intervals_as_intercept), ...
                'min_observations', 3);
            result = meg_fit_streaming_gls(metricDir, V_mask, blockTimes, ...
                double(cfg.frame_interval_seconds), cfg.contrasts, options);
            result.workflow = char(cfg.workflow);
            result.parent_workflow = char(cfg.parent_workflow);
            result.subject_id = subjectId;
            result.scan_id = scanId;
            result.parent_branch = char(cfg.parent_branch);
            result.parent_parameter_id = char(cfg.parent_parameter_id);
            result.parent_parameters = cfg.parent_parameters;
            result.block_table = blockFile;
            result.block_times = blockTimes;
            result.metric_directory = metricDir;
            result.metric_manifest_record = metricManifest;
            result.created = datestr(now,30);
            save(level1File, 'result', '-v7.3');
            manifest = level1_manifest(result, cfg);
            write_json(fullfile(level1Dir,'level1_scan_manifest.json'), manifest);
        end
        scanResults{iScan} = result;
        scanSummary(iScan).status = 'complete';
        scanSummary(iScan).frame_count = result.frame_count;
        scanSummary(iScan).condition_counts = result.condition_counts;
    catch ME
        scanSummary(iScan).status = 'failed';
        scanSummary(iScan).error_identifier = ME.identifier;
        scanSummary(iScan).error_message = getReport(ME,'extended','hyperlinks','off');
    end
end

participantFile = fullfile(subjectRoot, char(cfg.participant_subdir), ...
    'level1_participant.mat');
complete = strcmp({scanSummary.status}, 'complete');
if all(complete)
    if ~isequal(scanResults{1}.V_mask, scanResults{2}.V_mask) || ...
            ~isequal(scanResults{1}.mask_index, scanResults{2}.mask_index)
        error('MEG:WorkflowScanMask', 'Scan masks differ for %s.', subjectId);
    end
    ensure_directory(fileparts(participantFile));
    copeByScan = cat(4, scanResults{1}.cope, scanResults{2}.cope);
    participant = struct('schema_version',1, 'status','complete', ...
        'workflow',char(cfg.workflow), 'subject_id',subjectId, ...
        'scan_ids',{scanIds}, ...
        'scan_aggregation','equal_weight_mean_scan6_scan7', ...
        'pattern_names',{scanResults{1}.pattern_names}, ...
        'contrasts',cfg.contrasts, ...
        'contrast_vector',scanResults{1}.contrast_vector, ...
        'V_mask',scanResults{1}.V_mask, ...
        'mask_index',scanResults{1}.mask_index, ...
        'cope_by_scan',copeByScan, 'cope',mean(copeByScan,4), ...
        'parent_branch',char(cfg.parent_branch), ...
        'parent_parameter_id',char(cfg.parent_parameter_id), ...
        'created',datestr(now,30));
    save(participantFile, 'participant', '-v7.3');
    write_json(fullfile(fileparts(participantFile), ...
        'level1_participant_manifest.json'), struct('schema_version',1, ...
        'status','complete', 'workflow',char(cfg.workflow), ...
        'subject_id',subjectId, 'scan_ids',{scanIds}, ...
        'scan_aggregation',participant.scan_aggregation, ...
        'parent_branch',char(cfg.parent_branch), ...
        'parent_parameter_id',char(cfg.parent_parameter_id), ...
        'cope_dimensions',size(participant.cope), 'created',participant.created));
end

summary = struct('schema_version',1, 'status','complete', ...
    'workflow',char(cfg.workflow), 'subject_id',subjectId, ...
    'scan_ids',{scanIds}, 'scans',scanSummary, ...
    'participant_file',participantFile, 'config_file',configFile, ...
    'created',datestr(now,30));
if ~all(complete), summary.status = 'incomplete'; end
write_json(fullfile(subjectRoot,'subject_summary.json'), summary);
save(fullfile(subjectRoot,'subject_summary.mat'), 'summary');
if ~all(complete)
    error('MEG:WorkflowSubjectFailed', '%s failed: %d of 2 scans incomplete.', ...
        subjectId, nnz(~complete));
end
end

function cfg = load_config(filename)
if ~isfile(filename), error('MEG:WorkflowConfig','Missing config: %s',filename); end
cfg = jsondecode(fileread(filename));
required = {'workflow','parent_workflow','processed_root_env', ...
    'processed_root','results_root_env','results_root','scan_ids', ...
    'derived_subdir','cof_directory','metric_subdir','parent_results_subdir', ...
    'block_table_name','parent_branch','parent_parameter_id','parent_parameters', ...
    'results_subdir','level1_subdir','participant_subdir', ...
    'frame_interval_seconds','expected_blocks_per_condition_per_scan', ...
    'include_unlabeled_intervals_as_intercept','contrasts'};
for i=1:numel(required)
    if ~isfield(cfg,required{i}), error('MEG:WorkflowConfigField', ...
            'Missing config field %s.',required{i}); end
end
end

function validate_parent(cfg, resultsRoot, subjectId, scanId, metricDir, blockFile)
if ~isfolder(metricDir), error('MEG:WorkflowMetrics','Missing metrics: %s',metricDir); end
for name = {'divcurl_manifest.json','divcurl_statistics.mat'}
    if ~isfile(fullfile(metricDir,name{1})), error('MEG:WorkflowMetrics', ...
            'Missing %s.',fullfile(metricDir,name{1})); end
end
manifest = jsondecode(fileread(fullfile(metricDir,'divcurl_manifest.json')));
if ~isfield(manifest,'status') || ~strcmpi(char(manifest.status),'complete')
    error('MEG:WorkflowMetricsStatus','D/C metrics are incomplete: %s',metricDir);
end
lockFile = fullfile(resultsRoot,char(cfg.parent_results_subdir),subjectId, ...
    ['scan-' scanId],'cof_locked_manifest.json');
if ~isfile(lockFile), error('MEG:WorkflowParentLock','Missing %s.',lockFile); end
lock = jsondecode(fileread(lockFile));
if ~isfield(lock,'status') || ~strcmpi(char(lock.status),'complete') || ...
        ~strcmp(char(lock.workflow),char(cfg.parent_workflow)) || ...
        ~strcmp(char(lock.subject_id),subjectId) || ...
        ~strcmp(char(string(lock.scan_id)),scanId) || ...
        ~strcmp(char(lock.parameter_id),char(cfg.parent_parameter_id)) || ...
        ~any(strcmp(normalize_text_list(lock.branches),char(cfg.parent_branch)))
    error('MEG:WorkflowParentLock','Parent parameter/branch manifest mismatch: %s',lockFile);
end
if ~isfile(blockFile), error('MEG:WorkflowBlocks','Missing block table: %s',blockFile); end
end

function [V_mask, manifest] = load_metric_mask(metricDir)
loaded = load(fullfile(metricDir,'divcurl_statistics.mat'),'V_mask');
if ~isfield(loaded,'V_mask'), error('MEG:WorkflowMask','Missing V_mask.'); end
V_mask = logical(loaded.V_mask);
if ~isequal(size(V_mask),[200 500]) || ~any(V_mask(:))
    error('MEG:WorkflowMask','Expected nonempty 200x500 V_mask.');
end
manifest = jsondecode(fileread(fullfile(metricDir,'divcurl_manifest.json')));
end

function blocks = read_block_times(filename, subjectId, scanId, cfg)
T = readtable(filename);
required = {'subject_id','scan_id','branch','block_index','condition_code', ...
    'window_start_seconds','window_end_seconds'};
for i=1:numel(required)
    if ~ismember(required{i},T.Properties.VariableNames)
        error('MEG:WorkflowBlockField','%s lacks %s.',filename,required{i});
    end
end
if height(T) ~= 3*double(cfg.expected_blocks_per_condition_per_scan)
    error('MEG:WorkflowBlockCount','Expected 24 rows in %s.',filename);
end
if any(string(T.subject_id) ~= string(subjectId)) || ...
        any(string(T.scan_id) ~= string(scanId)) || ...
        any(string(T.branch) ~= string(cfg.parent_branch))
    error('MEG:WorkflowBlockIdentity','Block table identity mismatch: %s',filename);
end
codes = double(T.condition_code);
for code=1:3
    if nnz(codes==code) ~= double(cfg.expected_blocks_per_condition_per_scan)
        error('MEG:WorkflowBlockCount','Condition %d does not have 8 blocks.',code);
    end
end
starts = double(T.window_start_seconds); ends = double(T.window_end_seconds);
if any(~isfinite(starts)) || any(~isfinite(ends)) || any(ends<=starts)
    error('MEG:WorkflowBlockBounds','Invalid block bounds: %s',filename);
end
blocks = repmat(struct('block_index',0,'condition_code',0, ...
    'start_seconds',0,'end_seconds',0),1,height(T));
for i=1:height(T)
    blocks(i).block_index = double(T.block_index(i));
    blocks(i).condition_code = codes(i);
    blocks(i).start_seconds = starts(i);
    blocks(i).end_seconds = ends(i);
end
end

function validate_level1(result,V_mask,cfg,subjectId,scanId)
required={'cope','varcope','tstat','beta_gls','rho_hat','V_mask', ...
    'mask_index','pattern_names','contrasts','condition_counts', ...
    'workflow','subject_id','scan_id','parent_parameter_id'};
for i=1:numel(required)
    if ~isfield(result,required{i}), error('MEG:WorkflowLevel1Field', ...
            'Missing %s.',required{i}); end
end
if ~strcmp(char(result.workflow),char(cfg.workflow)) || ...
        ~strcmp(char(result.subject_id),subjectId) || ...
        ~strcmp(char(string(result.scan_id)),scanId) || ...
        ~strcmp(char(result.parent_parameter_id),char(cfg.parent_parameter_id)) || ...
        ~isequal(result.V_mask,V_mask) || ...
        ~isequal(size(result.cope),[nnz(V_mask),numel(cfg.contrasts),4])
    error('MEG:WorkflowLevel1Mismatch','Cached Level-1 result mismatch.');
end
end

function manifest = level1_manifest(result,cfg)
manifest=struct('schema_version',1,'status','complete', ...
    'workflow',char(cfg.workflow),'parent_workflow', ...
    char(cfg.parent_workflow),'subject_id',result.subject_id, ...
    'scan_id',result.scan_id,'solver',result.solver, ...
    'segment_aware',result.segment_aware,'parent_branch', ...
    char(cfg.parent_branch),'parent_parameter_id',char(cfg.parent_parameter_id), ...
    'parent_parameters',cfg.parent_parameters,'metric_directory', ...
    result.metric_directory,'block_table',result.block_table, ...
    'frame_interval_seconds',cfg.frame_interval_seconds, ...
    'pattern_names',{result.pattern_names},'contrasts',cfg.contrasts, ...
    'cope_dimensions',size(result.cope),'varcope_dimensions',size(result.varcope), ...
    'tstat_dimensions',size(result.tstat),'mask_size',size(result.V_mask), ...
    'vertex_count',nnz(result.V_mask),'frame_count',result.frame_count, ...
    'condition_counts',result.condition_counts,'created',result.created);
end

function tf=is_complete_level1(filename,cfg,subjectId,scanId)
tf=false; mf=fullfile(fileparts(filename),'level1_scan_manifest.json');
if ~isfile(filename)||~isfile(mf),return,end
try
    m=jsondecode(fileread(mf));
    tf=isfield(m,'status')&&strcmpi(char(m.status),'complete')&& ...
        strcmp(char(m.workflow),char(cfg.workflow))&& ...
        strcmp(char(m.subject_id),subjectId)&&strcmp(char(string(m.scan_id)),scanId)&& ...
        strcmp(char(m.parent_parameter_id),char(cfg.parent_parameter_id));
catch, tf=false; end
end

function root=resolve_root(cfg,envField,fallbackField)
root=getenv(char(cfg.(envField))); if isempty(root),root=char(cfg.(fallbackField));end
end
function ids=normalize_scan_ids(ids)
if isnumeric(ids),ids=arrayfun(@num2str,ids(:)','UniformOutput',false); ...
elseif ischar(ids),ids={ids}; elseif isstring(ids),ids=cellstr(ids(:)'); end
for i=1:numel(ids),ids{i}=validate_identifier(ids{i},'scanId');end
end
function id=validate_identifier(id,name)
id=char(string(id)); if isempty(regexp(id,'^[A-Za-z0-9_-]+$','once'))
    error('MEG:WorkflowIdentifier','Invalid %s.',name); end
end
function values=normalize_text_list(value)
if ischar(value),values={value};elseif isstring(value),values=cellstr(value(:)); ...
elseif iscell(value),values=cellfun(@char,value,'UniformOutput',false); ...
else,values=cellstr(string(value(:)));end
end
function ensure_directory(d),if ~isfolder(d),mkdir(d);end,end
function write_json(filename,value)
fid=fopen(filename,'w');if fid<0,error('MEG:WorkflowWrite','Cannot open %s.',filename);end
c=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value),'char');clear c
end
