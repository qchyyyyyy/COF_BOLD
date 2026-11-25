function rest_surrogate_velocity(i_batch, batch_size, cfg)
%CALC_VVF_HCP_YA_SURROGATE_FROM_PREPROCESSED Generate surrogate data and estimate VVFs.
%
% This script starts from existing TC_mesh_r_f_zscore.mat files. It does not
% rerun CIFTI preprocessing. For each subject, it writes:
%   <subject>/surrogate/TC_surrogate.mat
%   <subject>/surrogate/VelField_lh_COF.mat
%   <subject>/surrogate/VelField_rh_COF.mat
%
% Example:
%   rest_surrogate_velocity(1, 1)
%
% Optional configuration:
%   cfg.Type = 'rfMRI';
%   cfg.Session = {'REST1'};
%   cfg.PE = 'LR';
%   cfg.surrogateMethod = '2D'; % or '3D'

if nargin < 1 || isempty(i_batch)
    i_batch = 1;
end
if nargin < 2 || isempty(batch_size)
    batch_size = 1;
end
if nargin < 3 || isempty(cfg)
    cfg = struct();
end

repoRoot = find_repo_root(fileparts(mfilename('fullpath')));
cfg = set_default_cfg(cfg, repoRoot);

templateFile = fullfile(cfg.dataRoot, 'templates', 'coordinate_cii_v2.mat');
if ~exist(templateFile, 'file')
    error('Template file not found: %s', templateFile);
end
load(templateFile, 'cii_mask_left', 'cii_mask_right');

sessions = cfg.Session;
if ischar(sessions)
    sessions = {sessions};
end

for i_session = 1:numel(sessions)
    Experiment = [cfg.Type '_' sessions{i_session} '_' cfg.PE];
    save_dir = fullfile(cfg.dataRoot, cfg.Site, Experiment);

    if ~exist(save_dir, 'dir')
        warning('Preprocessed data directory not found: %s', save_dir);
        continue
    end

    subjectList = list_subjects_with_file(save_dir, cfg.inputFileName);
    if isempty(subjectList)
        warning('No subjects with %s found under %s', cfg.inputFileName, save_dir);
        continue
    end

    idx_sub = i_batch:batch_size:numel(subjectList);
    for i_sub = idx_sub
        subjectID = subjectList{i_sub};
        subject_dir = fullfile(save_dir, subjectID);
        surrogate_dir = fullfile(subject_dir, 'surrogate');
        input_file = fullfile(subject_dir, cfg.inputFileName);

        fprintf('Subject %d/%d: %s\n', i_sub, numel(subjectList), subjectID);
        fprintf('Generating %s surrogate data...\n', cfg.surrogateMethod);

        randomSeed = [];
        if ~isempty(cfg.randomSeed)
            randomSeed = cfg.randomSeed + i_sub;
        end

        is_valid = generate_surrogate(input_file, surrogate_dir, cfg.surrogateMethod, ...
            cii_mask_left, cii_mask_right, randomSeed);

        if ~is_valid
            continue
        end

        for i_method = 1:numel(cfg.Method)
            method = cfg.Method{i_method};
            fprintf('VVF estimation via %s...\n', method);

            vfCfg = [];
            vfCfg.path = ensure_trailing_filesep(surrogate_dir);
            vfCfg.fileName = cfg.surrogateFileName;

            fprintf('Left hemisphere...\n');
            vfCfg.hemisphere = 'lh';
            estimate_velocity(vfCfg, cii_mask_left, method, cfg.rho, cfg.tau, cfg.sigma);

            fprintf('Right hemisphere...\n');
            vfCfg.hemisphere = 'rh';
            estimate_velocity(vfCfg, cii_mask_right, method, cfg.rho, cfg.tau, cfg.sigma);
        end
    end
end

end

function cfg = set_default_cfg(cfg, repoRoot)

cfg = set_default(cfg, 'dataRoot', fullfile(repoRoot, 'data'));
cfg = set_default(cfg, 'Site', 'HCP-YA');
cfg = set_default(cfg, 'Type', 'rfMRI');
cfg = set_default(cfg, 'Session', {'REST1'});
cfg = set_default(cfg, 'PE', 'LR');
cfg = set_default(cfg, 'Method', {'COF'});
cfg = set_default(cfg, 'surrogateMethod', '3D');
cfg = set_default(cfg, 'inputFileName', 'TC_mesh_r_f_zscore.mat');
cfg = set_default(cfg, 'surrogateFileName', 'TC_surrogate.mat');
cfg = set_default(cfg, 'rho', 5e-4);
cfg = set_default(cfg, 'tau', 1e-8);
cfg = set_default(cfg, 'sigma', 8);
cfg = set_default(cfg, 'randomSeed', []);

end

function cfg = set_default(cfg, fieldName, defaultValue)

if ~isfield(cfg, fieldName) || isempty(cfg.(fieldName))
    cfg.(fieldName) = defaultValue;
end

end

function subjectList = list_subjects_with_file(save_dir, fileName)

items = dir(save_dir);
items = items([items.isdir]);
subjectList = {};

for i_item = 1:numel(items)
    name = items(i_item).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue
    end
    if exist(fullfile(save_dir, name, fileName), 'file')
        subjectList{end+1, 1} = name; %#ok<AGROW>
    end
end

end

function pathOut = ensure_trailing_filesep(pathIn)

pathOut = pathIn;
if isempty(pathOut) || pathOut(end) ~= filesep
    pathOut = [pathOut filesep];
end

end

function repoRoot = find_repo_root(startDir)

repoRoot = startDir;
while true
    if isfile(fullfile(repoRoot, 'setup.m'))
        return
    end
    parentDir = fileparts(repoRoot);
    if strcmp(parentDir, repoRoot)
        error('Could not locate repository root from: %s', startDir);
    end
    repoRoot = parentDir;
end

end
