function Data = behaviour_load_metrics(cfg)
%FIG3_BEHAVIORAL_LOAD_HCP_WM Load aligned HCP WM behaviour and D/C motif data.

groupDir = cfg.groupDir;
info = load(fullfile(groupDir, [cfg.experiment '_info.mat']), ...
    'Stim_Type_subs', 'Acc', 'RT');
SubjectList = behaviour_subject_list(fullfile(groupDir, [cfg.experiment '_SubList.mat']));

Acc = info.Acc;
RT = info.RT;
Stim_Type_subs = info.Stim_Type_subs;

if strcmp(cfg.session, 'WM')
    Acc(Acc == 0) = nan;
end

Subind = sum((squeeze(Stim_Type_subs(:,3,:)) - cfg.blockOrder).^2) == 0 & ...
    ~isnan(sum(Acc, 2))';
Stim_Type_subs = Stim_Type_subs(:,:,Subind);
SubjectList = SubjectList(Subind);

meanStimType = mean(Stim_Type_subs, 3);
meanStimType(:,2) = meanStimType(:,1) + meanStimType(:,2);

Events = zeros(cfg.eventLength, 1);
for iBlock = 1:size(Stim_Type_subs,1)
    idx = round(meanStimType(iBlock,1)):round(meanStimType(iBlock,2))-1;
    idx(idx > cfg.eventLength) = [];
    idx(idx < 1) = [];
    if numel(idx) > 38
        idx = idx(1:38);
    end
    Events(idx) = meanStimType(iBlock,3);
end

idx0Back = find(Events < 5 & Events > 0);
idx2Back = find(Events > 4);

AccRaw = Acc(Subind, cfg.blockColumns);
RTBlockValues = RT(Subind, cfg.blockColumns);
AccModelValues = AccRaw;
lambdaBoxCox = nan;
if cfg.useBoxCoxAccuracy
    assert(exist('boxcox','file')==2,'Behaviour:Dependency','Original boxcox function is required.');
    [accVec, lambdaBoxCox] = boxcox(AccModelValues(:));
    AccModelValues = reshape(accVec, [], numel(cfg.blockColumns));
end

nSub = nnz(Subind);
AccBlock = zeros(nSub, cfg.eventLength);
AccRawBlock = zeros(nSub, cfg.eventLength);
RTBlock = zeros(nSub, cfg.eventLength);
for iBlock = 1:numel(cfg.blockColumns)
    idxEvent = Events == iBlock;
    nEvent = nnz(idxEvent);
    AccBlock(:, idxEvent) = repmat(AccModelValues(:, iBlock), 1, nEvent);
    AccRawBlock(:, idxEvent) = repmat(AccRaw(:, iBlock), 1, nEvent);
    RTBlock(:, idxEvent) = repmat(RTBlockValues(:, iBlock), 1, nEvent);
end

subjectId = subjectListToNumeric(SubjectList);
behav = readtable(cfg.behaviourCsv);
[isFound, loc] = ismember(subjectId, behav.Subject);
if ~all(isFound)
    error('Could not align %d HCP subjects to behavioural CSV.', nnz(~isFound));
end
behav = behav(loc, :);

Age = parseHcpAge(behav.Age);
Gender01 = double(string(behav.Gender) == "M");

Num = loadGlobalFeaturePair(groupDir, cfg, 'num', Subind);
DCPQ = Num.Source + Num.Sink;

Size = loadGlobalFeaturePair(groupDir, cfg, 'size', Subind);
DCPS = (Size.Source + Size.Sink) ./ 2;

GlobalOP = loadOrderParamGlobal(groupDir, cfg, Subind);

Data = struct();
Data.SubjectList = SubjectList;
Data.Subject = subjectId;
Data.Subind = Subind;
Data.Events = Events;
Data.ConditionNames = cfg.conditionNames;
Data.ConditionIndex = {idx0Back, idx2Back};
Data.AccBlock = AccBlock;
Data.AccRawBlock = AccRawBlock;
Data.RTBlock = RTBlock;
Data.AccByBlock = AccModelValues;
Data.AccRawByBlock = AccRaw;
Data.RTByBlock = RTBlockValues;
Data.AccBoxCoxLambda = lambdaBoxCox;
Data.Age = Age;
Data.Gender01 = Gender01;
Data.BehaviorTable = behav;
Data.DCPQ = DCPQ;
Data.DCPS = DCPS;
Data.GlobalOP = GlobalOP;
end

function subjectId = subjectListToNumeric(SubjectList)
if isnumeric(SubjectList)
    subjectId = SubjectList(:);
    return
end
if isstring(SubjectList)
    subjectId = str2double(SubjectList(:));
    return
end
subjectId = str2double(string(SubjectList(:)));
end

function Age = parseHcpAge(ageColumn)
ageString = string(ageColumn);
Age = nan(numel(ageString), 1);
for iAge = 1:numel(ageString)
    parts = split(ageString(iAge), '-');
    vals = str2double(parts);
    vals = vals(isfinite(vals));
    if isempty(vals)
        Age(iAge) = 36;
    elseif numel(vals) == 1
        Age(iAge) = vals(1);
    else
        Age(iAge) = mean(vals);
    end
end
end

function P = loadGlobalFeaturePair(groupDir, cfg, featureName, Subind)
prefix = featurePrefix(featureName);
combinedFile = fullfile(groupDir, sprintf('%s_%sPattern_global_LR_%s.mat', ...
    cfg.experiment, featureName, cfg.method));
if exist(combinedFile, 'file')
    S = load(combinedFile);
    P.Source = S.([prefix 'Source'])(Subind,:);
    P.Sink = S.([prefix 'Sink'])(Subind,:);
    P.CCVortex = S.([prefix 'CCVortex'])(Subind,:);
    P.CVortex = S.([prefix 'CVortex'])(Subind,:);
    return
end

leftFile = fullfile(groupDir, sprintf('%s_%sPattern_global_%s_%s.mat', ...
    cfg.experiment, featureName, cfg.hemispheres{1}, cfg.method));
rightFile = fullfile(groupDir, sprintf('%s_%sPattern_global_%s_%s.mat', ...
    cfg.experiment, featureName, cfg.hemispheres{2}, cfg.method));
L = load(leftFile);
R = load(rightFile);
P.Source = L.([prefix 'Source'])(Subind,:) + R.([prefix 'Source'])(Subind,:);
P.Sink = L.([prefix 'Sink'])(Subind,:) + R.([prefix 'Sink'])(Subind,:);
P.CCVortex = L.([prefix 'CCVortex'])(Subind,:) + R.([prefix 'CCVortex'])(Subind,:);
P.CVortex = L.([prefix 'CVortex'])(Subind,:) + R.([prefix 'CVortex'])(Subind,:);
end

function prefix = featurePrefix(featureName)
switch lower(featureName)
    case 'num'
        prefix = 'Num';
    case 'size'
        prefix = 'Size';
    otherwise
        error('Unknown feature name: %s', featureName);
end
end

function GlobalOP = loadOrderParamGlobal(groupDir, cfg, Subind)
opFile = fullfile(groupDir, sprintf('%s_orderParam_global_LR_%s.mat', ...
    cfg.experiment, cfg.method));
S = load(opFile, 'order_param_global');
GlobalOP = S.order_param_global(Subind,:);
end
