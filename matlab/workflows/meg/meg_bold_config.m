function cfg = meg_bold_config(configFile)
%MEG_BOLD_CONFIG Load and validate the MEG-to-BOLD bridge contract.

if nargin ~= 1 || ~isfile(configFile)
    error('MEG:WorkflowConfig', 'Configuration file not found: %s', char(configFile));
end
cfg = jsondecode(fileread(configFile));
if ~isfield(cfg,'native_condition_subdir'), cfg.native_condition_subdir='meg_native_conditions'; end
required = {'workflow','processed_root_env','processed_root', ...
    'results_root_env','results_root','coordinate_file_env','coordinate_file', ...
    'source_subdirs','source_output_name_template','derived_subdir', ...
    'results_subdir','scan_ids','native_dt_seconds','output_tr_seconds', ...
    'hrf_parameters','maximum_imputed_hrf_mass','highpass_seconds', ...
    'wiener_lambda','lowpass_hz','source_block_size','mesh_time_block_size', ...
    'cof_time_block_frames','cof','metrics','inference'};
for i = 1:numel(required)
    if ~isfield(cfg, required{i})
        error('MEG:WorkflowConfigField', 'Missing configuration field %s.', required{i});
    end
end
if ~strcmp(char(cfg.workflow), 'meg_bold')
    error('MEG:WorkflowConfigID', 'Expected workflow meg_bold.');
end
cfg.processed_root = resolve_env(cfg, 'processed_root_env', 'processed_root');
cfg.results_root = resolve_env(cfg, 'results_root_env', 'results_root');
if isfield(cfg,'raw_root_env'), cfg.raw_root = resolve_env(cfg,'raw_root_env','raw_root'); end
if isfield(cfg,'fieldtrip_root_env'), cfg.fieldtrip_root = resolve_env(cfg,'fieldtrip_root_env','fieldtrip_root'); end
if isfield(cfg,'megconnectome_root_env'), cfg.megconnectome_root = resolve_env(cfg,'megconnectome_root_env','megconnectome_root'); end
if isfield(cfg,'spm_root_env'), cfg.spm_root = resolve_env(cfg,'spm_root_env','spm_root'); end
coordinateOverride = getenv(char(cfg.coordinate_file_env));
if ~isempty(coordinateOverride)
    cfg.coordinate_file = coordinateOverride;
elseif ~isfile(char(cfg.coordinate_file))
    codeRoot = fileparts(mfilename('fullpath'));
    error('NC:MEGCoordinate','Set COF_MEG_COORDINATE_FILE to the original MEG coordinate template.');
end
cfg.coordinate_file = char(cfg.coordinate_file);
if ~isfolder(cfg.processed_root)
    error('MEG:WorkflowProcessedRoot', 'Processed root not found: %s', cfg.processed_root);
end
if ~isfile(cfg.coordinate_file)
    error('MEG:WorkflowCoordinate', 'Coordinate file not found: %s', cfg.coordinate_file);
end
ratio = cfg.output_tr_seconds / cfg.native_dt_seconds;
if abs(ratio - round(ratio)) > 1e-10
    error('MEG:WorkflowSampling', 'output_tr_seconds must be an integer multiple of native_dt_seconds.');
end
if numel(cfg.hrf_parameters) ~= 7 || cfg.maximum_imputed_hrf_mass < 0 || ...
        cfg.maximum_imputed_hrf_mass > 1 || cfg.lowpass_hz <= 0
    error('MEG:WorkflowConfigValue', 'Invalid HRF, missingness, or filter configuration.');
end
if ~isfield(cfg.cof,'branches') || ~isfield(cfg.cof,'locked') || ...
        ~all(isfield(cfg.cof.locked,{'name','rho','tau','sigma'}))
    error('MEG:WorkflowCOFConfig','COF branches and named locked parameters are required.');
end
parameterName=char(cfg.cof.locked.name);
if isempty(regexp(parameterName,'^[A-Za-z0-9_-]+$','once'))
    error('MEG:WorkflowParameterName','Unsafe locked parameter name.');
end
if ~isfield(cfg.inference,'expected_subject_count') || ...
        ~isfield(cfg.inference,'required_scan_ids') || ...
        ~isfield(cfg.inference,'require_complete_cohort')
    error('MEG:WorkflowInferenceConfig','Complete-cohort inference fields are required.');
end
if ~isfolder(cfg.results_root), mkdir(cfg.results_root); end
end

function value = resolve_env(cfg, envField, valueField)
value = getenv(char(cfg.(envField)));
if isempty(value), value = char(cfg.(valueField)); end
end
