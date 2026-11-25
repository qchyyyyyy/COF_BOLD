function cfg = rest_config(overrides)
%FIG1_HCP3T_CONFIG Shared defaults for Fig. 1 HCP 3T fMRI validation scripts.
%
if nargin < 1 || isempty(overrides)
    overrides = struct();
end

repoRoot = repository_root();

cfg = struct();
cfg.repoRoot = repoRoot;
cfg.programRoot = repoRoot;
% Current phase: use the project-level data root. Keep this as the only
% switch point for the later reproduction-local data migration.
cfg.dataRoot = fullfile(cfg.programRoot, 'data');
cfg.processedDataRoot = fullfile(cfg.programRoot, 'data', 'processed');
cfg.resultsRoot = getenv('COF_OUTPUT_ROOT');
if isempty(cfg.resultsRoot), cfg.resultsRoot=fullfile(repoRoot,'results'); end
externalDataRoot=getenv('COF_DATA_ROOT');
if ~isempty(externalDataRoot), cfg.dataRoot=externalDataRoot; end
externalProcessedRoot=getenv('COF_PROCESSED_ROOT');
if ~isempty(externalProcessedRoot), cfg.processedDataRoot=externalProcessedRoot; end

cfg.Site = 'HCP-YA';
cfg.Type = 'rfMRI';
cfg.Session = 'REST1';
cfg.PE = 'LR';
cfg.Experiment = '';

cfg.Method = {'OF', 'PG', 'COF'};
cfg.modeMethods = {'COF'};
cfg.plotMethods = {'COF'};
cfg.primaryMethod = 'COF';
cfg.hemisphere = 'LR';
cfg.parcellation = 'Yeo2011_7networks';
cfg.modeParcellation = 'HCPex_360';

cfg.filtSigma = 3;
cfg.do_zscore = 1;
cfg.do_binarize = 0;
cfg.downSampRate = 3;

cfg.patternSourceNames = {'Source', 'Sink', 'CCVortex', 'CVortex'};
cfg.patternDisplayNames = {'diverging', 'converging', 'CC-spiral', 'C-spiral'};
cfg.methodColor = [189 215 238; 248 203 173; 197 224 180] ./ 255;

cfg.makeFigures = true;
cfg.saveFigures = true;
cfg.saveData = true;
cfg.closeFigures = false;
cfg.figureFormats = {'png', 'pdf'};

cfg.simSize = [200 200 100];
cfg.phaseSampleCount = 400;
cfg.amplitudeBins = 100;
cfg.waveOmega = 1;
cfg.waveNumber = 4;
cfg.streamTime = 10;
cfg.computeVelocityField = true;

cfg.modeMetric = 'All';
cfg.forceModeComputation = false;
cfg.makeGradientExampleFigure = true;
cfg.makeModeTimeFigure = true;
cfg.makeModeMapFigure = true;
cfg.exampleSubject = getenv('COF_EXAMPLE_SUBJECT');
cfg.exampleTimeWindow = 101:180;
cfg.modeTimeWindow = 101:180;
cfg.gradientFile = '';

cfg = apply_overrides(cfg, overrides);

if isempty(cfg.Experiment)
    sessionName = cfg.Session;
    if iscell(sessionName)
        sessionName = sessionName{1};
    end
    cfg.Experiment = [cfg.Type '_' sessionName '_' cfg.PE];
end

if contains(cfg.Type, 'rfMRI')
    defaultVvfDataDir = fullfile(cfg.dataRoot, 'processed', cfg.Site, cfg.Experiment);
else
    % Task-state velocity fields use the historical task tree directly;
    % resting-state derivatives are kept under data/processed.
    defaultVvfDataDir = fullfile(cfg.dataRoot, cfg.Site, cfg.Experiment);
end
cfg.vvfDataDir = set_default_path(cfg, overrides, 'vvfDataDir', defaultVvfDataDir);
cfg.participantMetricDir = set_default_path(cfg, overrides, 'participantMetricDir', ...
    fullfile(cfg.dataRoot, 'results', cfg.Site, cfg.Experiment, 'participant'));
cfg.groupMetricDir = set_default_path(cfg, overrides, 'groupMetricDir', ...
    fullfile(cfg.dataRoot, 'results', cfg.Site, cfg.Experiment, 'group'));
cfg.figureDir = set_default_path(cfg, overrides, 'figureDir', ...
    fullfile(cfg.resultsRoot, 'rest'));
cfg.figureDataDir = set_default_path(cfg, overrides, 'figureDataDir', ...
    fullfile(cfg.figureDir, 'data'));
cfg.gradientFile = set_default_path(cfg, overrides, 'gradientFile', ...
    fullfile(cfg.dataRoot, 'templates', 'gradients', 'g1.mat'));

end

function cfg = apply_overrides(cfg, overrides)

if ~isstruct(overrides)
    error('Overrides must be a struct.');
end

fields = fieldnames(overrides);
for i_field = 1:numel(fields)
    cfg.(fields{i_field}) = overrides.(fields{i_field});
end

end

function value = set_default_path(cfg, overrides, fieldName, defaultValue)

if isfield(overrides, fieldName) && ~isempty(overrides.(fieldName))
    value = overrides.(fieldName);
elseif isfield(cfg, fieldName) && ~isempty(cfg.(fieldName))
    value = cfg.(fieldName);
else
    value = defaultValue;
end

end
