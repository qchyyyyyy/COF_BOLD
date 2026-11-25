function Model = behaviour_metric_models(Data, cfg, opts)
%FIG3_BEHAVIORAL_FIT_METRIC_MODELS Fit appendix metric-behaviour models.
%
% Original age/sex-adjusted quadratic models for motif size and global OP.

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = fillMetricOptions(cfg, opts);

ensureDir(opts.OutputDir);

metrics = buildAppendixMetrics(Data, cfg, opts);
if isempty(metrics)
    error('No appendix metrics were available.');
end

outcomes = buildOutcomeByCondition(Data, cfg);
sourceRows = table();
statsRows = table();

for iMetric = 1:numel(metrics)
    statsStore = cell(numel(cfg.behaviors), numel(cfg.conditionNames));
    metric = metrics(iMetric);

    for iBeh = 1:numel(cfg.behaviors)
        behavior = cfg.behaviors{iBeh};
        Ycond = outcomes.(behavior);

        for iCond = 1:numel(cfg.conditionNames)
            condition = cfg.conditionNames{iCond};
            x = metric.ValuesByCondition(:, iCond);
            y = Ycond(:, iCond);
            st = fitQuadraticCovariateAdjusted(x, y, Data.Age, Data.Gender01, opts);
            st.Metric = string(metric.Name);
            st.MetricLabel = string(metric.Label);
            st.MetricCategory = string(metric.Category);
            st.Behavior = string(behavior);
            st.Condition = string(condition);
            statsStore{iBeh, iCond} = st;

            sourceRows = [sourceRows; table( ... %#ok<AGROW>
                repmat(string(metric.Name), numel(x), 1), ...
                repmat(string(metric.Label), numel(x), 1), ...
                repmat(string(metric.Category), numel(x), 1), ...
                repmat(string(condition), numel(x), 1), ...
                repmat(string(behavior), numel(x), 1), ...
                Data.Subject(:), x(:), y(:), Data.Age(:), Data.Gender01(:), ...
                st.FiniteCase(:), st.IncludedInModel(:), ...
                st.ExcludedAsInfluential(:), st.CooksDistance(:), ...
                'VariableNames', {'Metric', 'MetricLabel', 'MetricCategory', ...
                'Condition', 'Behavior', 'Subject', 'Predictor', 'Outcome', ...
                'Age', 'Gender01', 'FiniteCase', 'IncludedInModel', ...
                'ExcludedAsInfluential', 'CooksDistance'})];

            statsRows = [statsRows; packStatsRow(st, cfg)]; %#ok<AGROW>
        end
    end


end


sourceCsv = fullfile(opts.OutputDir, ...
    'behaviour_metric_source.csv');
statsCsv = fullfile(opts.OutputDir, ...
    'behaviour_metric_statistics.csv');
writetable(sourceRows, sourceCsv);
writetable(statsRows, statsCsv);

Model = struct();
Model.SourceCsv = sourceCsv;
Model.StatsCsv = statsCsv;
Model.Stats = statsRows;
Model.Source = sourceRows;
Model.Metrics = metrics;
Model.OutputDir = opts.OutputDir;
Model.Options = opts;
end

function opts = fillMetricOptions(cfg, opts)
opts.OutputDir = getOption(opts, 'OutputDir', cfg.appendixMetricDir);
opts.ModelMode = "ordinary";
opts.OutlierMethod = "none";
end

function value = getOption(opts, name, defaultValue)
if isfield(opts, name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end
end

function metrics = buildAppendixMetrics(Data, cfg, opts)
metrics = struct('Name', {}, 'Label', {}, 'Category', {}, 'ValuesByCondition', {});
metrics(end+1) = makeMetric('DCPS', 'D/C PS = (source + sink size)/2', ...
    'motif_size', conditionMeans(Data.DCPS, Data.ConditionIndex));
metrics(end+1) = makeMetric('GlobalOP', 'Global OP', ...
    'motif_order_parameter', conditionMeans(Data.GlobalOP, Data.ConditionIndex));

end

function metric = makeMetric(name, label, category, valuesByCondition)
metric = struct();
metric.Name = name;
metric.Label = label;
metric.Category = category;
metric.ValuesByCondition = valuesByCondition;
end

function values = conditionMeans(X, conditionIndex)
values = nan(size(X, 1), numel(conditionIndex));
for iCond = 1:numel(conditionIndex)
    values(:, iCond) = mean(X(:, conditionIndex{iCond}), 2, 'omitnan');
end
end

function outcomes = buildOutcomeByCondition(Data, cfg)
outcomes = struct();
for iBeh = 1:numel(cfg.behaviors)
    behavior = cfg.behaviors{iBeh};
    switch behavior
        case 'Acc'
            Ymat = Data.AccBlock;
        case 'RT'
            Ymat = Data.RTBlock;
        otherwise
            error('Unknown behavior: %s', behavior);
    end
    outcomes.(behavior) = conditionMeans(Ymat, Data.ConditionIndex);
end
end

function st = fitQuadraticCovariateAdjusted(x, y, age, gender01, opts)
finiteCase = isfinite(x) & isfinite(y) & isfinite(age) & isfinite(gender01);

st = struct();
st.NInitial = nnz(finiteCase);
st.NExcluded = 0;
st.N = st.NInitial;
st.ModelMode = opts.ModelMode;
st.OutlierMethod = opts.OutlierMethod;
st.OutlierThreshold = nan;
st.FiniteCase = finiteCase;
st.IncludedInModel = finiteCase;
st.ExcludedAsInfluential = false(size(x));
st.CooksDistance = nan(size(x));
st.IsValid = false;
st.XMean = nan;
st.XStd = nan;
st.XMin = nan;
st.XMax = nan;
st.OptimumMetric = nan;
st.BetaMetric = nan;
st.BetaMetric2 = nan;
st.P_Metric = nan;
st.P_Metric2 = nan;
st.RSquared = nan;
st.AdjustedRSquared = nan;
st.XGrid = [];
st.YFit = [];
st.YCI = [];
st.X = x(finiteCase);
st.Y = y(finiteCase);

if st.NInitial < 10 || std(x(finiteCase), 'omitnan') <= 0 || ...
        std(y(finiteCase), 'omitnan') <= 0
    return
end

xFit = x(st.IncludedInModel);
yFitInput = y(st.IncludedInModel);
ageFit = age(st.IncludedInModel);
genderFit = gender01(st.IncludedInModel);

st.N = numel(xFit);
st.NExcluded = st.NInitial - st.N;
st.X = xFit;
st.Y = yFitInput;
st.IsValid = st.N >= 10 && std(xFit, 'omitnan') > 0 && ...
    std(yFitInput, 'omitnan') > 0;

if ~st.IsValid
    return
end

[mdl, xMean, xStd] = fitQuadraticCore(xFit, yFitInput, ageFit, genderFit);
st.XMean = xMean;
st.XStd = xStd;
st.XMin = min(xFit);
st.XMax = max(xFit);

coefNames = string(mdl.CoefficientNames);
idxMetric = coefNames == "MetricZ";
idxMetric2 = coefNames == "MetricZ2";
st.BetaMetric = mdl.Coefficients.Estimate(idxMetric);
st.BetaMetric2 = mdl.Coefficients.Estimate(idxMetric2);
st.P_Metric = mdl.Coefficients.pValue(idxMetric);
st.P_Metric2 = mdl.Coefficients.pValue(idxMetric2);
st.RSquared = mdl.Rsquared.Ordinary;
st.AdjustedRSquared = mdl.Rsquared.Adjusted;
st.Model = mdl;

if isfinite(st.BetaMetric2) && abs(st.BetaMetric2) > eps
    optimumZ = -st.BetaMetric ./ (2 * st.BetaMetric2);
    st.OptimumMetric = st.XMean + optimumZ .* st.XStd;
end

pad = max(1e-6, rangeOrOne(xFit) * 0.04);
xGrid = linspace(min(xFit) - pad, max(xFit) + pad, 160)';
zGrid = (xGrid - st.XMean) ./ st.XStd;
pred = table( ...
    zeros(numel(xGrid), 1), zGrid, zGrid.^2, ...
    repmat(mean(ageFit, 'omitnan'), numel(xGrid), 1), ...
    repmat(mean(genderFit, 'omitnan'), numel(xGrid), 1), ...
    'VariableNames', {'Y', 'MetricZ', 'MetricZ2', 'Age', 'Gender'});
[yFit, yCI] = predict(mdl, pred, 'Prediction', 'curve');
st.XGrid = xGrid;
st.YFit = yFit;
st.YCI = yCI;
end

function [mdl, xMean, xStd] = fitQuadraticCore(x, y, age, gender01)
xMean = mean(x, 'omitnan');
xStd = std(x, 'omitnan');
MetricZ = (x - xMean) ./ xStd;
MetricZ2 = MetricZ.^2;
Age = age(:);
Gender = gender01(:);
Y = y(:);
T = table(Y, MetricZ, MetricZ2, Age, Gender);
mdl = fitlm(T, 'Y ~ MetricZ + MetricZ2 + Age + Gender');
end

function row = packStatsRow(st, cfg)
formula = sprintf('Outcome ~ z(%s) + z(%s)^2 + Age + Gender', ...
    st.MetricLabel, st.MetricLabel);
row = table( ...
    st.Metric, st.MetricLabel, st.MetricCategory, st.Condition, st.Behavior, ...
    st.NInitial, st.NExcluded, st.N, st.ModelMode, st.OutlierMethod, ...
    st.OutlierThreshold, st.BetaMetric, st.BetaMetric2, st.P_Metric, st.P_Metric2, ...
    st.RSquared, st.AdjustedRSquared, st.OptimumMetric, string(formula), ...
    'VariableNames', {'Metric', 'MetricLabel', 'MetricCategory', 'Condition', ...
    'Behavior', 'NInitial', 'NExcluded', 'N', 'ModelMode', 'OutlierMethod', ...
    'OutlierThreshold', 'Beta_MetricZ', 'Beta_MetricZ2', 'P_Metric', ...
    'P_Metric2', 'RSquared', 'AdjustedRSquared', 'Optimum_Metric', ...
    'ModelFormula'});
if st.Behavior == "Acc"
    row.OutcomeLabel = string(cfg.behaviorLabels('Acc'));
else
    row.OutcomeLabel = string(cfg.behaviorLabels('RT'));
end
end

function r = rangeOrOne(x)
r = max(x) - min(x);
if ~isfinite(r) || r <= 0
    r = 1;
end
end

function ensureDir(pathValue)
if ~exist(pathValue, 'dir')
    mkdir(pathValue);
end
end
