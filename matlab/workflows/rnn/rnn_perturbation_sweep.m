function R = rnn_perturbation_sweep(cfg)
% Compute the six manuscript perturbation families for each model and trial.
% PQ is mean(source + sink count); PS is mean(source + sink size).
% The order of retained design rows preserves seeds of the six-family config.

%% Defaults

if nargin < 1 || isempty(cfg)
    cfg = struct();
end

assert(isfield(cfg,'rootDir') && isfolder(cfg.rootDir),'RNN:InputRoot','Set cfg.rootDir to the trained models.');
cfg = setDefault(cfg, 'matFilePattern', 'activation_dynamics_baseline_1296_%s.mat');
cfg = setDefault(cfg, 'subjects', arrayfun(@(i) sprintf('sub%02d', i), 1:100, 'UniformOutput', false));
assert(isfield(cfg,'outDir') && ~isempty(cfg.outDir),'RNN:OutputRoot','Set a separate output directory.');
assert(~strcmp(cfg.rootDir,cfg.outDir),'RNN:OutputRoot','Input and output directories must differ.');
cfg = setDefault(cfg, 'trialIndices', 1:3);
cfg = setDefault(cfg, 'force', true);

cfg = setDefault(cfg, 'xFields', {'x_features_2back','x_features','X_features','X','x'});
cfg = setDefault(cfg, 'yFields', {'true_targets_2back','true_targets','targets','Y','y'});
cfg = setDefault(cfg, 'h0Fields', {'h0_2back','h0','h0_used'});

cfg = setDefault(cfg, 'inputNoiseDims', 'stimulus_only');   % 'all', 'stimulus_only', numeric vector
cfg = setDefault(cfg, 'randomSeed', 20260610);
cfg = setDefault(cfg, 'nRepeats', 1);
cfg = setDefault(cfg, 'saveHiddenSamples', false);

% Severity grids. You can narrow these for pilot runs.

% Stress-like perturbations
cfg = setDefault(cfg, 'inputNoiseSigmas',  [0.15 0.25 0.35 0.45 0.60]);
cfg = setDefault(cfg, 'inputDropoutRates', [0.05 0.10 0.20 0.30 0.45]);
cfg = setDefault(cfg, 'hiddenNoiseSigmas', [0.10 0.20 0.40 0.80]);

% Recurrent structural perturbations
cfg = setDefault(cfg, 'recurrentPruneFractions', [0.25 0.30 0.35 0.40 0.45]);
cfg = setDefault(cfg, 'recurrentNoiseStds',      [0.20 0.40 0.70 1.00]);

% Leak / time-scale perturbations
cfg = setDefault(cfg, 'leakGlobalScales',    [0.75 0.65 0.55 0.45 0.35]);

% Gain / excitability perturbations

% Low-rank recurrent perturbation

if ~exist(cfg.outDir, 'dir')
    mkdir(cfg.outDir);
end

if exist('rnn_hidden_metrics', 'file') ~= 2
    error(['Cannot find rnn_hidden_metrics.m on MATLAB path. ', ...
           'Please add your project functions folder first.']);
end

%% Build perturbation design
Design = buildPerturbationDesign(cfg);
fprintf('Perturbation design contains %d method/severity/repeat rows.\n', height(Design));

%% Main loop
AllRows = table();
HiddenSamples = struct();

for iSub = 1:numel(cfg.subjects)
    subID = cfg.subjects{iSub};
    matFile = resolveMatFile(cfg, subID);
    if ~exist(matFile, 'file')
        error('RNN:MissingInput', 'Missing model file for %s: %s.', subID, matFile);
    end

    fprintf('\n[%d/%d] Loading %s\n', iSub, numel(cfg.subjects), matFile);
    S = load(matFile);

    try
        params0 = loadRNNParamsFromStruct(S);
        Xall = getFirstExisting(S, cfg.xFields, 'stimulus features');
        Yall = getFirstExisting(S, cfg.yFields, 'labels/targets');
        H0all = getFirstExisting(S, cfg.h0Fields, 'initial hidden state');
    catch ME
        rethrow(ME);
    end

    coords = [];
    if isfield(S, 'coordinates_2d_FD')
        coords = S.coordinates_2d_FD;
    elseif isfield(S, 'xy_coords')
        coords = S.xy_coords;
    else
        error('RNN:MissingCoordinates', 'No coordinates for %s.', subID);
    end

    trialList = resolveTrialList(cfg.trialIndices, Xall);
    for it = 1:numel(trialList)
        trialIdx = trialList(it);
        fprintf('  Trial %d\n', trialIdx);

        try
            xClean = extractTrialStimulus(Xall, trialIdx);
            yTrue  = extractTrialLabels(Yall, trialIdx, size(xClean,1));
            h0     = extractTrialH0(H0all, trialIdx, params0.N);
        catch ME
            rethrow(ME);
        end

        [Hclean, logitsClean] = rolloutLeakyRNNMatlabExtended(xClean, h0, params0, struct());
        behClean = computeBehaviorMetrics(logitsClean, yTrue);
        patClean = safeHiddenToPatternMetrics(Hclean, coords);
        rowClean = makeResultRow(subID, iSub, trialIdx, 'clean', 'none', 0, 0, 0, ...
            behClean, patClean, behClean, patClean);
        AllRows = [AllRows; rowClean]; %#ok<AGROW>

        for iD = 1:height(Design)
            d = Design(iD,:);
            seed = cfg.randomSeed + 1000000*iSub + 1000*trialIdx + iD;
            rng(seed, 'twister');

            try
                params = params0;
                xPert = xClean;
                rolloutOpt = struct();

                switch string(d.method)
                    case "input_noise"
                        xPert = addStimulusNoise(xClean, d.severity, cfg);
                    case "input_dropout"
                        xPert = applyInputDropout(xClean, d.severity, cfg);
                    case "hidden_noise"
                        rolloutOpt.hiddenNoiseSigma = d.severity;
                        rolloutOpt.hiddenNoiseSeed = seed;
                    case "recurrent_prune"
                        params.W_hh = applyRecurrentPrune(params0.W_hh, d.severity);
                    case "recurrent_noise"
                        params.W_hh = applyRecurrentNoise(params0.W_hh, d.severity);
                    case "leak_global"
                        params.alpha_vec = clampAlpha(params0.alpha_vec .* d.severity);
                    otherwise
                        error('Unknown perturbation method: %s', string(d.method));
                end

                [Hpert, logitsPert] = rolloutLeakyRNNMatlabExtended(xPert, h0, params, rolloutOpt);
                behPert = computeBehaviorMetrics(logitsPert, yTrue);
                patPert = safeHiddenToPatternMetrics(Hpert, coords);

                row = makeResultRow(subID, iSub, trialIdx, string(d.method), string(d.family), ...
                    d.severity, d.fraction, d.repeat, behPert, patPert, behClean, patClean);
                AllRows = [AllRows; row]; %#ok<AGROW>

                if cfg.saveHiddenSamples
                    key = matlab.lang.makeValidName(sprintf('%s_trial%d_%s_%g_rep%d', subID, trialIdx, string(d.method), d.severity, d.repeat));
                    HiddenSamples.(key).H = single(Hpert);
                    HiddenSamples.(key).logits = single(logitsPert);
                end

            catch ME
                rethrow(ME);
            end
        end
    end
end

%% Save
csvFile = fullfile(cfg.outDir, 'rnn_perturbation_sweep_all_methods.csv');
matFile = fullfile(cfg.outDir, 'rnn_perturbation_sweep_all_methods.mat');
writetable(AllRows, csvFile);
save(matFile, 'AllRows', 'Design', 'cfg', 'HiddenSamples', '-v7.3');

fprintf('\nSaved table: %s\n', csvFile);
fprintf('Saved MAT:   %s\n', matFile);


R = struct();
R.AllRows = AllRows;
R.Design = Design;
R.cfg = cfg;
R.csvFile = csvFile;
R.matFile = matFile;

end

%% Design helpers
function D = buildPerturbationDesign(cfg)
methods = {};
families = {};
severities = [];
fractions = [];
repeats = [];

for rep = 1:cfg.nRepeats
    for v = cfg.inputNoiseSigmas(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'input_noise', 'stress', v, NaN, rep);
    end
    for v = cfg.inputDropoutRates(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'input_dropout', 'stress', v, NaN, rep);
    end
    for v = cfg.hiddenNoiseSigmas(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'hidden_noise', 'stress', v, NaN, rep);
    end
    for v = cfg.recurrentPruneFractions(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'recurrent_prune', 'structural', v, v, rep);
    end
    for v = cfg.recurrentNoiseStds(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'recurrent_noise', 'structural', v, NaN, rep);
    end
    for v = cfg.leakGlobalScales(:)'
        [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, 'leak_global', 'structural', v, NaN, rep);
    end
end

D = table(string(methods), string(families), severities, fractions, repeats, ...
    'VariableNames', {'method','family','severity','fraction','repeat'});
end

function [methods, families, severities, fractions, repeats] = addDesignRow(methods, families, severities, fractions, repeats, method, family, sev, frac, rep)
methods{end+1,1} = method;
families{end+1,1} = family;
severities(end+1,1) = sev;
fractions(end+1,1) = frac;
repeats(end+1,1) = rep;
end

%% RNN loading / rollout
function params = loadRNNParamsFromStruct(S)
W_hh = getFirstExisting(S, {'W_hh','Whh','W_rec'}, 'recurrent weight W_hh');
N = size(W_hh,1);
if size(W_hh,2) ~= N
    error('W_hh must be square. Got %s.', mat2str(size(W_hh)));
end

W_ih = getFirstExisting(S, {'W_ih','W_in','W_input'}, 'input weight W_ih');
if size(W_ih,1) ~= N && size(W_ih,2) == N
    W_ih = W_ih';
end
if size(W_ih,1) ~= N
    error('W_ih shape %s is incompatible with N=%d.', mat2str(size(W_ih)), N);
end

b = getFirstExisting(S, {'b','bias','b_h'}, 'hidden bias b');
b = b(:);
if numel(b) ~= N
    error('hidden bias length %d incompatible with N=%d.', numel(b), N);
end

W_out = getFirstExisting(S, {'W_out','fc_out_weight','W_readout'}, 'readout weight W_out');
if size(W_out,2) ~= N && size(W_out,1) == N
    W_out = W_out';
end
if size(W_out,2) ~= N
    error('W_out shape %s is incompatible with N=%d.', mat2str(size(W_out)), N);
end

fc_out_bias = getFirstExistingOptional(S, {'fc_out_bias','b_out','out_bias','fc_bias'});
if isempty(fc_out_bias)
    fc_out_bias = 0;
end
fc_out_bias = double(fc_out_bias(:));
if isempty(fc_out_bias); fc_out_bias = 0; end

alpha = getFirstExistingOptional(S, {'alpha','leak_alpha'});
if isempty(alpha)
    alpha = 0.2;
end
alpha = double(alpha);
if isscalar(alpha)
    alpha_vec = repmat(alpha, N, 1);
else
    alpha_vec = alpha(:);
    if numel(alpha_vec) ~= N
        error('alpha has %d entries; expected scalar or N=%d.', numel(alpha_vec), N);
    end
end

params = struct();
params.W_hh = double(W_hh);
params.W_ih = double(W_ih);
params.b = double(b);
params.W_out = double(W_out);
params.fc_out_bias = double(fc_out_bias(1));
params.alpha = double(mean(alpha_vec));
params.alpha_vec = clampAlpha(alpha_vec);
params.gain_vec = ones(N,1);
params.N = N;
params.F = size(W_ih,2);
end

function [H, logits] = rolloutLeakyRNNMatlabExtended(X, h0, params, opt)
if nargin < 4; opt = struct(); end
opt = setDefault(opt, 'hiddenNoiseSigma', 0);
opt = setDefault(opt, 'hiddenNoiseSeed', []);

X = double(X);
[T,F] = size(X);
N = params.N;
if F ~= params.F
    error('Input feature dimension mismatch: X has F=%d, W_ih expects F=%d.', F, params.F);
end

h = double(h0(:));
if numel(h) ~= N
    error('h0 has %d entries; expected N=%d.', numel(h), N);
end

if ~isempty(opt.hiddenNoiseSeed)
    rng(opt.hiddenNoiseSeed, 'twister');
end

alpha_vec = params.alpha_vec(:);
gain_vec = params.gain_vec(:);
H = zeros(N,T);
logits = zeros(T,1);

for t = 1:T
    pre = params.W_ih * X(t,:)' + params.W_hh * h + params.b;
    pre = gain_vec .* pre;
    h = (1 - alpha_vec) .* h + alpha_vec .* tanh(pre);
    if opt.hiddenNoiseSigma > 0
        h = h + opt.hiddenNoiseSigma * randn(N,1);
    end
    H(:,t) = h;
    logits(t,1) = params.W_out * h + params.fc_out_bias;
end
end

%% Perturbations
function Xn = addStimulusNoise(X, sigma, cfg)
Xn = double(X);
if sigma <= 0
    return;
end
dims = resolveInputNoiseDims(cfg.inputNoiseDims, size(X,2));
Xn(:,dims) = Xn(:,dims) + sigma * randn(size(Xn,1), numel(dims));
end

function Xd = applyInputDropout(X, rate, cfg)
Xd = double(X);
dims = resolveInputNoiseDims(cfg.inputNoiseDims, size(X,2));
mask = rand(size(Xd,1), numel(dims)) > rate;
Xd(:,dims) = Xd(:,dims) .* mask;
end

function dims = resolveInputNoiseDims(setting, F)
if isnumeric(setting)
    dims = setting(:)';
elseif ischar(setting) || isstring(setting)
    s = lower(string(setting));
    if s == "all"
        dims = 1:F;
    elseif s == "stimulus_only"
        if F <= 2
            dims = 1:F;
        else
            dims = 1:(F-2);
        end
    else
        error('Unknown inputNoiseDims setting: %s', s);
    end
else
    error('inputNoiseDims must be numeric, char, or string.');
end
dims = dims(dims >= 1 & dims <= F);
if isempty(dims)
    error('No valid input noise dimensions resolved.');
end
end

function W = applyRecurrentPrune(W0, fraction)
W = double(W0);
if fraction <= 0; return; end
mask = rand(size(W)) > fraction;
W = W .* mask;
end

function W = applyRecurrentNoise(W0, noiseStd)
W0 = double(W0);
baseStd = std(W0(:), 'omitnan');
if ~isfinite(baseStd) || baseStd <= 0
    baseStd = 1;
end
W = W0 + noiseStd * baseStd * randn(size(W0));
end

function alpha_vec = clampAlpha(alpha_vec)
alpha_vec = double(alpha_vec(:));
alpha_vec = min(max(alpha_vec, 0.001), 0.999);
end

%% Behavior and pattern metrics
function metrics = computeBehaviorMetrics(logits, yTrue)
z = double(logits(:));
y = double(yTrue(:));
assert(numel(z)==numel(y),'RNN:Length','Predictions and labels have different lengths.');
p = 1 ./ (1 + exp(-z));
pred = p >= 0.5;
ybin = y >= 0.5;
eps0 = 1e-8;
bce = -mean(ybin .* log(p + eps0) + (1-ybin) .* log(1-p + eps0), 'omitnan');
signed = 2*double(ybin) - 1;
margin = mean(signed .* z, 'omitnan');

metrics = struct();
metrics.accuracy = mean(pred == ybin, 'omitnan');
metrics.bce = bce;
metrics.margin = margin;
metrics.hit_rate = mean(pred(ybin == 1) == 1, 'omitnan');
metrics.false_alarm_rate = mean(pred(ybin == 0) == 1, 'omitnan');
metrics.target_acc = mean(pred(ybin == 1) == ybin(ybin == 1), 'omitnan');
metrics.nontarget_acc = mean(pred(ybin == 0) == ybin(ybin == 0), 'omitnan');
end

function pat = safeHiddenToPatternMetrics(H, coords)
try
    if isempty(coords)
        error('Missing coordinates_2d_FD.');
    end
    raw = rnn_hidden_metrics(H, coords);
    pat = summarizePatternMetrics(raw);
catch ME
    rethrow(ME);
end
end

function pat = summarizePatternMetrics(raw)
names = {'NumSource','NumSink','SizeSource','SizeSink','order_param_global'};
pat = struct();
for i = 1:numel(names)
    nm = names{i};
    if isfield(raw, nm)
        arr = double(raw.(nm)(:));
    else
        error('RNN:MetricField','Missing required pattern metric: %s.',nm);
    end
    pat.([nm '_mean']) = mean(arr, 'omitnan');
    pat.([nm '_std']) = std(arr, 'omitnan');
    pat.([nm '_sum']) = sum(arr, 'omitnan');
    pat.([nm '_median']) = median(arr, 'omitnan');
end

pat.PQ_source_mean = pat.NumSource_mean;
pat.PQ_sink_mean = pat.NumSink_mean;
pat.PQ_DC_mean = pat.NumSource_mean + pat.NumSink_mean;
pat.PS_source_mean = pat.SizeSource_mean;
pat.PS_sink_mean = pat.SizeSink_mean;
pat.PS_DC_mean = pat.SizeSource_mean + pat.SizeSink_mean;
end

function row = makeResultRow(subID, subNum, trialIdx, method, family, severity, fraction, repeat, beh, pat, behClean, patClean)
row = table();
row.subject = string(subID);
row.subject_num = subNum;
row.trial_index = trialIdx;
row.method = string(method);
row.family = string(family);
row.severity = severity;
row.fraction = fraction;
row.repeat = repeat;

row.accuracy = beh.accuracy;
row.bce = beh.bce;
row.margin = beh.margin;
row.hit_rate = beh.hit_rate;
row.false_alarm_rate = beh.false_alarm_rate;
row.target_acc = beh.target_acc;
row.nontarget_acc = beh.nontarget_acc;

row.clean_accuracy = behClean.accuracy;
row.clean_bce = behClean.bce;
row.clean_margin = behClean.margin;
row.delta_accuracy = beh.accuracy - behClean.accuracy;
row.delta_bce = beh.bce - behClean.bce;
row.delta_margin = beh.margin - behClean.margin;

patFields = fieldnames(pat);
for i = 1:numel(patFields)
    f = patFields{i};
    row.(f) = pat.(f);
    if isfield(patClean, f)
        row.(['clean_' f]) = patClean.(f);
        row.(['delta_' f]) = pat.(f) - patClean.(f);
    end
end
end

%% Plotting
function matFile = resolveMatFile(cfg, subID)
candidate1 = fullfile(cfg.rootDir, subID, sprintf(cfg.matFilePattern, subID));
candidate2 = fullfile(cfg.rootDir, sprintf(cfg.matFilePattern, subID));
candidate3 = sprintf(cfg.matFilePattern, subID);
if exist(candidate1, 'file')
    matFile = candidate1;
elseif exist(candidate2, 'file')
    matFile = candidate2;
else
    matFile = candidate3;
end
end

function val = getFirstExisting(S, candidates, label)
val = [];
for i = 1:numel(candidates)
    key = candidates{i};
    if isfield(S, key)
        val = S.(key);
        return;
    end
end
error('Could not find %s. Tried: %s', label, strjoin(candidates, ', '));
end

function val = getFirstExistingOptional(S, candidates)
val = [];
for i = 1:numel(candidates)
    key = candidates{i};
    if isfield(S, key)
        val = S.(key);
        return;
    end
end
end

function trials = resolveTrialList(setting, X)
if ischar(setting) || isstring(setting)
    if any(strcmpi(string(setting), ["all","*"]))
        trials = 1:inferNTrials(X);
    else
        error('RNN:Trials','Use a numeric trial-index array or the string all.');
    end
else
    trials = setting;
end
trials = trials(:)';
end

function n = inferNTrials(X)
X = squeeze(X);
if ndims(X) == 3
    s = size(X);
    if s(1) <= 256 && s(2) > s(1)
        n = s(1);
    elseif s(3) <= 256
        n = s(3);
    else
        n = 1;
    end
else
    n = 1;
end
end

function Xtrial = extractTrialStimulus(Xall, trialIdx)
X = squeeze(Xall);
if ndims(X) == 2
    assert(trialIdx==1,'RNN:Trials','A two-dimensional stimulus contains one trial.');
    Xtrial = double(X);
elseif ndims(X) == 3
    s = size(X);
    if s(1) <= 256 && s(2) > s(1)
        Xtrial = squeeze(X(trialIdx,:,:));       % [Trials,T,F]
    elseif s(3) <= 256
        Xtrial = squeeze(X(:,:,trialIdx));       % [T,F,Trials]
    else
        Xtrial = squeeze(X(:,:,trialIdx));
    end
else
    error('Unsupported X dimensions: %s', mat2str(size(Xall)));
end

if size(Xtrial,1) <= 32 && size(Xtrial,2) > 32
    Xtrial = Xtrial';
end
end

function y = extractTrialLabels(Yall, trialIdx, T)
Y = squeeze(Yall);
if isvector(Y)
    y = double(Y(:));
elseif ndims(Y) == 2
    s = size(Y);
    if any(s == T)
        if s(1) == T
            y = double(Y(:, trialIdx));
        else
            y = double(Y(trialIdx, :)');
        end
    else
        y = double(Y(:));
    end
elseif ndims(Y) == 3
    s = size(Y);
    if s(1) <= 256 && s(2) == T
        y = squeeze(Y(trialIdx,:,1));
    elseif s(1) == T
        y = squeeze(Y(:,1,trialIdx));
    else
        y = squeeze(Y(trialIdx,:,:));
    end
    y = double(y(:));
else
    error('Unsupported Y dimensions: %s', mat2str(size(Yall)));
end

assert(numel(y)==T,'RNN:Length','Labels and stimulus have different lengths.');
end

function h0 = extractTrialH0(H0all, trialIdx, N)
H0 = squeeze(H0all);
if isvector(H0)
    h0 = double(H0(:));
elseif ndims(H0) == 2
    if size(H0,1) == N
        if size(H0,2) == 1
            h0 = double(H0(:,1));
        else
            h0 = double(H0(:,trialIdx));
        end
    elseif size(H0,2) == N
        h0 = double(H0(trialIdx,:)');
    else
        error('H0 shape %s incompatible with N=%d.', mat2str(size(H0)), N);
    end
elseif ndims(H0) == 3
    if size(H0,1) == N
        h0 = double(H0(:,1,trialIdx));
    elseif size(H0,2) == N
        h0 = double(H0(trialIdx,:,1)');
    else
        error('H0 shape %s incompatible with N=%d.', mat2str(size(H0)), N);
    end
else
    error('Unsupported H0 dimensions: %s', mat2str(size(H0all)));
end

if numel(h0) ~= N
    error('Extracted h0 has %d entries, expected N=%d.', numel(h0), N);
end
end

function cfg = setDefault(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end
