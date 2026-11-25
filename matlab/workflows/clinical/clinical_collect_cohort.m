function D=clinical_collect_cohort(cfg,dataset)
% Source-faithful cohort reconstruction; known historical row exclusions retained.
assert(all(isfield(cfg,{'dataRoot','resultsRoot','method','hemisphere','excludedRows'})),'Missing paths/settings');
switch dataset
case 'ZIB',D=collectZIBData(cfg);
case 'ds000115',D=collectDS000115Data(cfg);
otherwise,error('Unknown dataset');
end
end

function D = collectZIBData(cfg)
    Site = 'SCZ_ZIB';
    data_path = fullfile(cfg.dataRoot, Site);
    results_path = fullfile(cfg.resultsRoot,Site);

    subjFile = fullfile(results_path, 'subjectList.mat');
    if isfield(cfg,'subjectListFile'),subjFile=cfg.subjectListFile;end
    if ~exist(subjFile,'file')
        error('Cannot find %s', subjFile);
    end
    load(subjFile, 'subjectList');

    % Historical row exclusion is explicit in the configuration.
    bad_id = cfg.excludedRows.ZIB;
    assert(all(bad_id>=1 & bad_id<=numel(subjectList)),'clinical:Exclusion','Invalid ZIB exclusion rows.');
    allSubjectList=string(subjectList(:));
    subjectList(bad_id) = [];

    groupLabel = strings(numel(subjectList),1);
    groupLabel(contains(subjectList,'HC')) = "HC";
    groupLabel(contains(subjectList,'SCZ')) = "SCZ";
    assert(all(ismember(groupLabel,["HC","SCZ"])),'clinical:Group','Unrecognised participant group.');
    group01 = double(groupLabel == "SCZ");

    % Demographics: age and sex are used only for group-effect annotation.
    [Age, Sex] = loadZIBDemographics(data_path, subjectList);
    covMat = [Age(:), Sex(:)];
    covNames = {'Age','Sex'};

    eventFile = fullfile(data_path, 'nback_Events_group.mat');
    if ~exist(eventFile,'file')
        error('Cannot find %s', eventFile);
    end
    load(eventFile, 'Events', 'Acc'); %#ok<LOAD>
    Events(bad_id,:) = [];
    Acc(bad_id,:) = [];

    patternFile = fullfile(results_path, 'group', ...
        sprintf('task-nback_dir-PA_run-1_numPattern_global_%s_%s.mat', cfg.hemisphere, cfg.method));
    if ~exist(patternFile,'file')
        error('Cannot find %s', patternFile);
    end
    P = load(patternFile);
    verify_pattern_order(P,allSubjectList,bad_id,cfg);
    NumSource = P.NumSource;
    NumSink = P.NumSink;
    assert(isequal(size(NumSource),size(NumSink)),'clinical:Shape','Source/sink matrices disagree.');
    NumSource(bad_id,:) = [];
    NumSink(bad_id,:) = [];
    DCPQ = NumSource + NumSink;

    condNames = ["0-back", "2-back"];
    condCodes = [1, 2];
    PQ = nan(numel(subjectList), numel(condNames));
    ACC = nan(numel(subjectList), numel(condNames));

    for c = 1:numel(condNames)
        M = Events == condCodes(c);
        PQ(:,c) = maskedMeanRows(DCPQ, M);
        ACC(:,c) = maskedMeanRows(Acc, M);
    end

    D = struct();
    D.dataset = 'ZIB';
    D.prettyName = 'ZIB';
    assert(numel(unique(string(subjectList)))==numel(subjectList),'clinical:IDs','Duplicate participant IDs.');
    D.subjectList = string(subjectList(:));
    D.groupLabel = groupLabel;
    D.group01 = group01(:);
    D.conditionNames = condNames;
    D.PQ = PQ;
    D.Acc = ACC;
    D.covMat = covMat;
    D.covNames = covNames;
end

function D = collectDS000115Data(cfg)
    Site = 'SCZ_WM';
    results_path = fullfile(cfg.resultsRoot,Site);

    partFile = fullfile(results_path, 'participants.tsv');
    if ~exist(partFile,'file')
        error('Cannot find %s', partFile);
    end
    subject_info = readtable(partFile, 'FileType','text', 'Delimiter','\t');
    numericFields={'age','yrschool','nback0_targ','nback0_nont','nback1_targ','nback1_nont','nback2_targ','nback2_nont'};
    for k=1:numel(numericFields)
        key=numericFields{k};
        if ~isnumeric(subject_info.(key))
            subject_info.(key)=str2double(string(subject_info.(key)));
        end
    end

    condit = string(subject_info.condit);
    idx = condit == "SCZ" | condit == "CON";
    subject_info = subject_info(idx,:);

    % Exclusion is applied after selecting SCZ/CON, as in the original code.
    bad_id=cfg.excludedRows.ds000115;
    assert(all(bad_id>=1 & bad_id<=height(subject_info)),'clinical:Exclusion','Invalid ds000115 exclusion rows.');
    allSubjectList=string(subject_info.participant_id);
    subject_info(bad_id,:) = [];

    subjectList = string(subject_info.participant_id);
    groupLabel = strings(height(subject_info),1);
    groupLabel(string(subject_info.condit) == "CON") = "HC";
    groupLabel(string(subject_info.condit) == "SCZ") = "SCZ";
    assert(all(ismember(groupLabel,["HC","SCZ"])),'clinical:Group','Unrecognised participant group.');
    group01 = double(groupLabel == "SCZ");

    Age = subject_info.age;
    Edu = subject_info.yrschool;
    sexText=upper(string(subject_info.gender));
    assert(all(ismember(sexText,["MALE","FEMALE"])),'clinical:Sex','Unknown sex coding.');
    Sex = double(sexText == "MALE");
    covMat = [Age(:), Edu(:), Sex(:)];
    covNames = {'Age','Edu','Sex'};

    patternFile = fullfile(results_path, 'group', ...
        sprintf('numPattern_global_%s_%s.mat', cfg.hemisphere, cfg.method));
    if ~exist(patternFile,'file')
        error('Cannot find %s', patternFile);
    end
    P = load(patternFile);
    verify_pattern_order(P,allSubjectList,bad_id,cfg);
    NumSource = P.NumSource;
    NumSink = P.NumSink;
    assert(isequal(size(NumSource),size(NumSink)),'clinical:Shape','Source/sink matrices disagree.');

    % Robust alignment in case the pattern matrix still contains the removed row.
    if size(NumSource,1) == height(subject_info) + numel(bad_id)
        NumSource(bad_id,:) = [];
        NumSink(bad_id,:) = [];
    elseif size(NumSource,1) ~= height(subject_info)
        error('NC:ClinicalAlignment','Pattern and participant counts disagree; positional truncation is forbidden.');
    end

    DCPQ = NumSource + NumSink;

    % Reconstruct event labels as in your previous ds000115 script.
    T = 475;
    Events = zeros(1,T);
    Events([38:38+146, 250:250+146]) = 1;
    Events = [Events, Events*2, Events*3];

    condNames = ["0-back", "1-back", "2-back"];
    condCodes = [1, 2, 3];
    PQ = nan(height(subject_info), numel(condNames));
    ACC = nan(height(subject_info), numel(condNames));

    Acc_0b = mean([subject_info.nback0_targ, subject_info.nback0_nont],2,'omitnan');
    Acc_1b = mean([subject_info.nback1_targ, subject_info.nback1_nont],2,'omitnan');
    Acc_2b = mean([subject_info.nback2_targ, subject_info.nback2_nont],2,'omitnan');
    ACC(:,1) = Acc_0b;
    ACC(:,2) = Acc_1b;
    ACC(:,3) = Acc_2b;

    for c = 1:numel(condNames)
        PQ(:,c) = maskedMeanRows(DCPQ, Events == condCodes(c));
    end

    D = struct();
    D.dataset = 'ds000115';
    D.prettyName = 'ds000115';
    assert(numel(unique(subjectList))==numel(subjectList),'clinical:IDs','Duplicate participant IDs.');
    D.subjectList = subjectList;
    D.groupLabel = groupLabel;
    D.group01 = group01(:);
    D.conditionNames = condNames;
    D.PQ = PQ;
    D.Acc = ACC;
    D.covMat = covMat;
    D.covNames = covNames;
end

function [Age, Sex] = loadZIBDemographics(data_path, subjectList)
    Age = nan(numel(subjectList),1);
    Sex = nan(numel(subjectList),1);

    hcFile = fullfile(data_path, 'HC_Demographics.csv');
    sczFile = fullfile(data_path, 'Patients_Demographics.csv');
    if ~exist(hcFile,'file') || ~exist(sczFile,'file')
        error('clinical:Demographics','Both ZIB demographic files are required.');
    end

    HC_Demographics = readtable(hcFile);
    SCZ_Demographics = readtable(sczFile);

    ids = cellfun(@(x) x(5:end), cellstr(subjectList), 'UniformOutput', false);
    hcKeep = ismember(string(HC_Demographics.SubjectID), string(ids));
    sczKeep = ismember(string(SCZ_Demographics.SubID), string(ids));
    HC_Demographics = HC_Demographics(hcKeep,:);
    SCZ_Demographics = SCZ_Demographics(sczKeep,:);

    demoSubjects = [strcat("sub-", string(HC_Demographics.SubjectID)); ...
                    strcat("sub-", string(SCZ_Demographics.SubID))];
    demoAge = [HC_Demographics.Age; SCZ_Demographics.Age];
    demoSex = [HC_Demographics.Sex; SCZ_Demographics.Sex];
    demoSex = parseSexToBinary(demoSex);

    assert(numel(unique(demoSubjects))==numel(demoSubjects),'clinical:IDs','Duplicate demographic IDs.');
    [tf, loc] = ismember(string(subjectList), demoSubjects);
    assert(all(tf),'clinical:IDs','Missing demographic participant.');
    Age(tf) = demoAge(loc(tf));
    Sex(tf) = demoSex(loc(tf));
end

function sex01 = parseSexToBinary(x)
    if isnumeric(x)
        assert(numel(unique(x(isfinite(x))))<=2,'clinical:Sex','Expected at most two original numeric sex codes.');
        sex01 = double(x);
        return;
    end
    xs = upper(string(x));
    assert(all(ismissing(xs)|ismember(xs,["M","MALE","1","F","FEMALE","0"])),'clinical:Sex','Unknown textual sex coding.');
    sex01 = nan(numel(xs),1);
    sex01(xs == "M" | xs == "MALE" | xs == "1") = 1;
    sex01(xs == "F" | xs == "FEMALE" | xs == "0") = 0;
end

function out = maskedMeanRows(X, M)
    X = double(X);
    if isvector(M)
        M = M(:)';
        M = repmat(M, size(X,1), 1);
    end
    M = logical(M);
    if ~isequal(size(X), size(M))
        error('Mask size [%s] does not match data size [%s].', ...
            num2str(size(M)), num2str(size(X)));
    end
    Xmasked = X;
    Xmasked(~M) = nan;
    out = mean(Xmasked, 2, 'omitnan');
end

function verify_pattern_order(P,original,excluded,cfg)
    retained=original;retained(excluded)=[];
    if size(P.NumSource,1)==numel(original),expected=original;
    elseif size(P.NumSource,1)==numel(retained),expected=retained;
    else,error('clinical:Alignment','Pattern row count differs from both original and retained cohorts.');end
    if isfield(P,'subjectList')
        assert(isequal(string(P.subjectList(:)),expected(:)),'clinical:Alignment','Pattern IDs/order differ from the cohort; positional relabelling is forbidden.');
    else
        assert(isfield(cfg,'allow_historical_positional_input') && cfg.allow_historical_positional_input,'clinical:Alignment','No pattern IDs; explicitly enable the documented original positional input contract.');
    end
end
