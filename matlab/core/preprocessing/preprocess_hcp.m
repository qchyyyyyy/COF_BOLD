function outputs = preprocess_hcp(config_file)
%CALC_VVF_HCP_YA Reproduce HCP-YA preprocessing and COF velocity fields.
%
% This repository entry point is derived from:
%   historical COF_BOLD/preprocess_hcp.m
%
% The numerical preprocessing and COF algorithms are preserved, while data,
% dependency, subject, and output paths are supplied by a JSON configuration.

pipeline_root = fileparts(mfilename('fullpath'));
repo_root = repository_root();
assert(nargin==1 && ~isempty(config_file),'HCP:Config','Supply an explicit preprocessing JSON configuration.');
config_file = char(config_file);
cfg = jsondecode(fileread(config_file));
validate_config(cfg, config_file);

addpath(cfg.spm_root);

resource_root = getenv('COF_TEMPLATE_ROOT');
if isfield(cfg,'template_root') && ~isempty(cfg.template_root), resource_root=char(cfg.template_root); end
assert(~isempty(resource_root),'Set cfg.template_root or COF_TEMPLATE_ROOT to the external template directory.');
mask_file = fullfile(resource_root, 'coordinate_cii_v2.mat');
surface_root = fullfile(resource_root, 'HCP_S1200_Atlas_Z4_pkXDZ');
load(mask_file, 'cii_mask_left', 'cii_mask_right');

surf_left = gifti(fullfile(surface_root, ...
    'S1200.L.flat.32k_fs_LR.surf.gii'));
surf_right = gifti(fullfile(surface_root, ...
    'S1200.R.flat.32k_fs_LR.surf.gii'));
x_left = double(surf_left.vertices(:, 1));
y_left = double(surf_left.vertices(:, 2));
x_right = double(surf_right.vertices(:, 1));
y_right = double(surf_right.vertices(:, 2));

x_target_cood = cfg.grid.x_min:cfg.grid.spacing:cfg.grid.x_max;
y_target_cood = cfg.grid.y_min:cfg.grid.spacing:cfg.grid.y_max;
subjects = normalize_text_list(cfg.subjects);
methods = normalize_text_list(cfg.methods);
run_specs = build_run_specs(cfg, subjects);
outputs = cell(numel(run_specs), 1);
output_index = 0;
for i_run = 1:numel(run_specs)
        subject = run_specs(i_run).subject;
        experiment = run_specs(i_run).experiment;
        input_file = run_specs(i_run).input_file;
        output_dir = run_specs(i_run).output_dir;

        fprintf('subject=%s experiment=%s\n', subject, experiment);
        fprintf('input=%s\n', input_file);
        fprintf('output=%s\n', output_dir);
        if ~exist(input_file, 'file')
            error('preprocess_hcp:MissingInput', ...
                'Required CIFTI input does not exist: %s', input_file);
        end
        if exist(output_dir, 'dir') && ~cfg.allow_existing_output
            existing = dir(fullfile(output_dir, '*.mat'));
            if ~isempty(existing)
                error('preprocess_hcp:ExistingOutput', ...
                    ['Output directory already contains MAT files. Refusing to ' ...
                     'overwrite without allow_existing_output=true: %s'], output_dir);
            end
        end

        if strcmpi(run_specs(i_run).mode, 'rfMRI')
            is_valid = preprocess_rest(input_file, output_dir, cfg.tr, ...
                cii_mask_left, cii_mask_right, x_target_cood, y_target_cood, ...
                x_left, y_left, x_right, y_right);
        else
            is_valid = preprocess_task(input_file, output_dir, cfg.tr, ...
                cii_mask_left, cii_mask_right, x_target_cood, y_target_cood, ...
                x_left, y_left, x_right, y_right);
        end
        if ~is_valid
            error('preprocess_hcp:PreprocessingFailed', ...
                'Preprocessing failed for %s %s.', subject, experiment);
        end

        for i_method = 1:numel(methods)
            method = methods{i_method};
            fprintf('method=%s hemisphere=lh\n', method);
            vvf_cfg = struct();
            vvf_cfg.path = output_dir;
            vvf_cfg.fileName = 'TC_mesh_r_f_zscore.mat';
            vvf_cfg.hemisphere = 'lh';
            estimate_velocity(vvf_cfg, cii_mask_left, method, ...
                cfg.cof.rho, cfg.cof.tau, cfg.cof.sigma);

            fprintf('method=%s hemisphere=rh\n', method);
            vvf_cfg.hemisphere = 'rh';
            estimate_velocity(vvf_cfg, cii_mask_right, method, ...
                cfg.cof.rho, cfg.cof.tau, cfg.cof.sigma);
        end

        output_index = output_index + 1;
        outputs{output_index} = output_dir;
end
outputs = outputs(1:output_index);
end

function values = normalize_text_list(value)
if ischar(value)
    values = {value};
elseif isstring(value)
    values = cellstr(value);
elseif iscell(value)
    values = cellfun(@char, value, 'UniformOutput', false);
else
    error('preprocess_hcp:InvalidConfig', ...
        'Expected a character, string, or cell array of text values.');
end
end

function validate_config(cfg, config_file)
required = {'output_root', 'spm_root', 'subjects', ...
    'phase_encoding', 'tr', 'methods', 'grid', 'cof', ...
    'allow_existing_output'};
for i_field = 1:numel(required)
    if ~isfield(cfg, required{i_field})
        error('preprocess_hcp:InvalidConfig', ...
            'Missing field "%s" in %s.', required{i_field}, config_file);
    end
end
if ~isfield(cfg, 'mode')
    cfg.mode = 'tfMRI';
end
if strcmpi(cfg.mode, 'rfMRI')
    if ~isfield(cfg, 'raw_rest_root')
        error('preprocess_hcp:InvalidConfig', ...
            'rfMRI configuration requires raw_rest_root.');
    elseif ~exist(cfg.raw_rest_root, 'dir')
        error('preprocess_hcp:InvalidConfig', ...
            'raw_rest_root does not exist: %s', cfg.raw_rest_root);
    end
    if ~isfield(cfg, 'sessions') || isempty(cfg.sessions)
        error('preprocess_hcp:InvalidConfig', ...
            'rfMRI configuration requires a nonempty sessions field.');
    end
elseif strcmpi(cfg.mode, 'tfMRI')
    if ~isfield(cfg, 'raw_task_root')
        error('preprocess_hcp:InvalidConfig', ...
            'tfMRI configuration requires raw_task_root.');
    elseif ~exist(cfg.raw_task_root, 'dir')
        error('preprocess_hcp:InvalidConfig', ...
            'raw_task_root does not exist: %s', cfg.raw_task_root);
    end
    if ~isfield(cfg, 'tasks') || isempty(cfg.tasks)
        error('preprocess_hcp:InvalidConfig', ...
            'tfMRI configuration requires a nonempty tasks field.');
    end
else
    error('preprocess_hcp:InvalidConfig', ...
        'mode must be tfMRI or rfMRI.');
end
if ~exist(cfg.spm_root, 'dir')
    error('preprocess_hcp:InvalidConfig', ...
        'spm_root does not exist: %s', cfg.spm_root);
end
if cfg.grid.spacing <= 0
    error('preprocess_hcp:InvalidConfig', 'grid.spacing must be positive.');
end
if cfg.tr <= 0
    error('preprocess_hcp:InvalidConfig', 'tr must be positive.');
end
if cfg.cof.sigma < 0
    error('preprocess_hcp:InvalidConfig', 'cof.sigma must be nonnegative.');
end
end

function specs = build_run_specs(cfg, subjects)
if ~isfield(cfg, 'mode')
    mode = 'tfMRI';
else
    mode = char(cfg.mode);
end
specs = struct('subject', {}, 'experiment', {}, 'input_file', {}, ...
    'output_dir', {}, 'mode', {});
index = 0;
if strcmpi(mode, 'rfMRI')
    sessions = normalize_text_list(cfg.sessions);
    for i_session = 1:numel(sessions)
        experiment = sprintf('rfMRI_%s_%s', sessions{i_session}, cfg.phase_encoding);
        for i_subject = 1:numel(subjects)
            subject = subjects{i_subject};
            index = index + 1;
            specs(index).subject = subject;
            specs(index).experiment = experiment;
            specs(index).input_file = fullfile(cfg.raw_rest_root, subject, ...
                'MNINonLinear', 'Results', experiment, ...
                [experiment '_Atlas_hp2000_clean.dtseries.nii']);
            specs(index).output_dir = fullfile(cfg.output_root, experiment, subject);
            specs(index).mode = 'rfMRI';
        end
    end
else
    tasks = cfg.tasks;
    for i_task = 1:numel(tasks)
        task = tasks(i_task);
        experiment = sprintf('tfMRI_%s_%s', task.session, cfg.phase_encoding);
        task_root = fullfile(cfg.raw_task_root, task.source_name, 'preproc');
        for i_subject = 1:numel(subjects)
            subject = subjects{i_subject};
            index = index + 1;
            specs(index).subject = subject;
            specs(index).experiment = experiment;
            specs(index).input_file = fullfile(task_root, subject, ...
                'MNINonLinear', 'Results', experiment, ...
                [experiment '_Atlas_MSMAll.dtseries.nii']);
            specs(index).output_dir = fullfile(cfg.output_root, experiment, subject);
            specs(index).mode = 'tfMRI';
        end
    end
end
end
