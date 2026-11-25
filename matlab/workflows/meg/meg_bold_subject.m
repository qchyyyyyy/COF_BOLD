function summary = meg_bold_subject(configFile, subjectId, scanId, stage)
%MEG_BOLD_SUBJECT MEG-power to BOLD-like signal, COF, and PQ/PS.
%
% stage is signal or full; full includes the manuscript COF parameters.

if nargin < 3 || nargin > 4
    error('MEG:WorkflowArguments', ...
        'Usage: meg_bold_subject(config, subject_id, scan_id, [stage])');
end
if nargin < 4 || isempty(stage), stage = 'full'; end
stage = char(stage);
if ~any(strcmp(stage,{'signal','full'}))
    error('MEG:WorkflowStage', 'stage must be signal or full.');
end
subjectId = validate_identifier(subjectId,'subject_id');
scanId = validate_identifier(scanId,'scan_id');
cfg = meg_bold_config(configFile);
codeRoot = fileparts(mfilename('fullpath'));
addpath(codeRoot);

sourceDir = resolve_source_dir(cfg,subjectId,scanId);
outputDir = fullfile(cfg.processed_root,subjectId,'Wrkmem', ...
    char(cfg.derived_subdir),['scan-' scanId]);
if ~isfolder(outputDir), mkdir(outputDir); end
resultDir = fullfile(cfg.results_root,char(cfg.results_subdir),subjectId,['scan-' scanId]);
if ~isfolder(resultDir), mkdir(resultDir); end

preprocessedFile = fullfile(outputDir,'preprocessed_signal.mat');
generatedFile = fullfile(outputDir,'generated_bold_source.mat');
if ~isfile(preprocessedFile)
    [power,validMask,~,sourceManifest] = meg_load_power_session(sourceDir);
    if abs(double(sourceManifest.frame_interval_seconds)-cfg.native_dt_seconds) > 1e-10
        error('MEG:WorkflowInputTR', 'Source frame interval differs from configuration.');
    end
    generated = meg_bold_generate(power,validMask,cfg);
    temporal = meg_bold_preprocess_source(generated,cfg);
    coordinates = load(cfg.coordinate_file,'x_left','y_left','x_right','y_right', ...
        'cii_mask_left','cii_mask_right');
    nLeft = numel(coordinates.x_left); nRight = numel(coordinates.x_right);
    if nLeft+nRight ~= size(power,1)
        error('MEG:WorkflowCoordinateCount', 'Coordinate/source count mismatch.');
    end
    exactL = meg_interpolate_hemisphere(temporal.exact_project1_source(1:nLeft,:), ...
        coordinates.x_left,coordinates.y_left,coordinates.cii_mask_left);
    exactR = meg_interpolate_hemisphere(temporal.exact_project1_source(nLeft+1:end,:), ...
        coordinates.x_right,coordinates.y_right,coordinates.cii_mask_right);
    data = struct;
    data.BrainImg_lh_exact_project1 = meg_grid_zscore(exactL,coordinates.cii_mask_left);
    data.BrainImg_rh_exact_project1 = meg_grid_zscore(exactR,coordinates.cii_mask_right);
    data.F_TR = cfg.output_tr_seconds;
    data.F_band = [1/cfg.highpass_seconds cfg.lowpass_hz];
    data.Normalization = 'gridwise temporal zscore';
    data.R_globalSig = 'yes';
    data.R_deconv_exact_project1 = 'yes';
    data.reliable_frame_mask = logical(generated.reliable_frame_mask);
    data.imputed_hrf_mass_fraction = generated.imputed_hrf_mass_fraction;
    data.time_seconds = generated.output_time_seconds;
    save(preprocessedFile,'data','-v7.3');
    save(generatedFile,'generated','temporal','-v7.3');
    clear power exactL exactR temporal generated
else
    loaded = load(preprocessedFile,'data'); data = loaded.data;
end

signalManifest = struct;
signalManifest.schema_version = 1;
signalManifest.status = 'complete';
signalManifest.workflow = 'meg_bold';
signalManifest.analysis_type = 'exploratory';
signalManifest.subject_id = subjectId;
signalManifest.scan_id = scanId;
signalManifest.stage = 'signal';
signalManifest.source_directory = sourceDir;
signalManifest.preprocessed_file = preprocessedFile;
signalManifest.branches = {'exact_project1'};
signalManifest.signal_branches = {'exact_project1'};
signalManifest.cof_branches = cfg.cof.branches;
signalManifest.output_frames = numel(data.time_seconds);
signalManifest.reliable_frames = nnz(data.reliable_frame_mask);
signalManifest.maximum_imputed_hrf_mass = cfg.maximum_imputed_hrf_mass;
signalManifest.hrf_parameters = cfg.hrf_parameters;
signalManifest.created = datestr(now,30);
write_json(fullfile(resultDir,'signal_manifest.json'),signalManifest);

if strcmp(stage,'signal')
    summary = signalManifest;
    return
end

summary = meg_bold_velocity(configFile,subjectId,scanId,'locked');
end

function sourceDir = resolve_source_dir(cfg,subjectId,scanId)
name = char(cfg.source_output_name_template);
name = strrep(name,'{subject_id}',subjectId);
name = strrep(name,'{scan_id}',scanId);
subdirs = cfg.source_subdirs;
if ischar(subdirs), subdirs={subdirs}; elseif isstring(subdirs), subdirs=cellstr(subdirs); end
for i=1:numel(subdirs)
    candidate = fullfile(cfg.processed_root,subjectId,'Wrkmem',char(subdirs{i}),name);
    if isfile(fullfile(candidate,'session_manifest.json'))
        sourceDir = candidate; return
    end
end
error('MEG:WorkflowSourceMissing','No compatible source-power directory for %s scan %s.',subjectId,scanId);
end

function identifier = validate_identifier(identifier,name)
identifier = char(string(identifier));
if isempty(regexp(identifier,'^[A-Za-z0-9_-]+$','once'))
    error('MEG:WorkflowIdentifier','%s contains unsupported characters.',name);
end
end

function write_json(filename,value)
fid=fopen(filename,'w');
if fid<0,error('MEG:WorkflowJSONWrite','Cannot open %s.',filename);end
cleanup=onCleanup(@()fclose(fid)); fwrite(fid,jsonencode(value),'char'); clear cleanup
end
