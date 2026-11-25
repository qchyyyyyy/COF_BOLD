function outputs=behaviour_group_profiles(cfg)
if ischar(cfg) || isstring(cfg), cfg=jsondecode(fileread(cfg)); end
if ~isfield(cfg,'verifyHistorical'), cfg.verifyHistorical=true; end
if cfg.verifyHistorical
    behaviour_verify_inputs(cfg.groupDir,fullfile(repository_root(),'configs','behaviour','historical_input_hashes.json'));
end
% Historical exploratory estimator; no plots, no overwrite of source caches.
assert(all(isfield(cfg,{'groupDir','behaviorCsv','outputDir'})),'Explicit inputs/output required');
assert(~strcmp(cfg.groupDir,cfg.outputDir),'Outputs must be separate from inputs');
if ~isfolder(cfg.outputDir),mkdir(cfg.outputDir);end
hemisphere={'lh','rh'};Method={'OF','PG','COF'};i_method=3;
Session='WM';Experiment='tfMRI_WM_LR';
results_path=[char(cfg.groupDir) filesep];behav_csv_path=char(cfg.behaviorCsv);
% Arm definition around fitted D/C optimum.
ARM_DELTA_Z = 0.25;

core_trait_specs = get_core_plot_specs();

%% ---------------- Load task timing and behavior ----------------
SubjectList = behaviour_subject_list([results_path Experiment '_SubList.mat']);
load([results_path Experiment '_info.mat'],'Stim_Type_subs','Acc','RT');

if strcmp(Session,'WM')
    Acc(Acc==0) = nan;
end

Subind = sum((squeeze(Stim_Type_subs(:,3,:))-[8;1;6;4;5;7;2;3]).^2)==0 & ~isnan(sum(Acc,2))';
Stim_Type_subs = Stim_Type_subs(:,:,Subind);
meanStimType = mean(Stim_Type_subs,3);
meanStimType(:,2) = meanStimType(:,1) + meanStimType(:,2);

Delay = 0;
T = 404;
Events = zeros(T,1);
for i_block = 1:size(Stim_Type_subs,1)
    idx = round(meanStimType(i_block,1)):round(meanStimType(i_block,2))-1;
    idx = idx + Delay;
    idx(idx>404)=[]; idx(idx<1)=[];
    if numel(idx) > 38
        idx = idx(1:38);
    end
    Events(idx) = meanStimType(i_block,3);
end
idx_0B = find(Events<5 & Events>0);
idx_2B = find(Events>4);

idx_block = 1:8;
Acc_raw = Acc(Subind,idx_block);
assert(exist('boxcox','file')==2,'Behaviour:Dependency','Original boxcox function is required.');
[Acc_boxcox_vec, lambda] = boxcox(Acc_raw(:)); %#ok<ASGLU>
Acc_boxcox = reshape(Acc_boxcox_vec,[],length(idx_block));
RT = RT(Subind,idx_block);

Acc_block_boxcox = nan([nnz(Subind) T]);
RT_block = Acc_block_boxcox;
for i_block = 1:8
    for i_sub = 1:nnz(Subind)
        Acc_block_boxcox(i_sub,Events==i_block) = Acc_boxcox(i_sub,i_block);
        RT_block(i_sub,Events==i_block) = RT(i_sub,i_block);
    end
end

SubjectList = SubjectList(Subind);
subj_num = str2double(cellstr(SubjectList));

%% ---------------- Load and align HCP unrestricted behavioral traits ----------------
behav_info_all = readtable(behav_csv_path);
[is_match, loc] = ismember(subj_num, behav_info_all.Subject);
assert(all(is_match), 'Behaviour:MissingParticipant', 'Behaviour CSV lacks an imaging participant.');

keep = is_match;
SubjectList = SubjectList(keep);
subj_num = subj_num(keep);
Acc_raw = Acc_raw(keep,:);
Acc_boxcox = Acc_boxcox(keep,:);
Acc_block_boxcox = Acc_block_boxcox(keep,:);
RT = RT(keep,:);
RT_block = RT_block(keep,:);

behav_info = behav_info_all(loc(keep),:);
sub_gender = behav_info.Gender;
sub_gender_bi = zeros(height(behav_info),1);
sub_gender_bi(strcmp(sub_gender,'M')) = 1;
sub_age = parse_hcp_age(behav_info.Age);

%% ---------------- Compute D/C PQ and optimum ----------------
Num_L = load([results_path Experiment '_numPattern_global_' hemisphere{1} '_' Method{i_method} '.mat']);
Num_R = load([results_path Experiment '_numPattern_global_' hemisphere{2} '_' Method{i_method} '.mat']);

NumSource = Num_L.NumSource(Subind,:) + Num_R.NumSource(Subind,:);
NumSink   = Num_L.NumSink(Subind,:)   + Num_R.NumSink(Subind,:);
NumSource = NumSource(keep,:);
NumSink   = NumSink(keep,:);

PQ_source_2b = mean(NumSource(:,idx_2B),2,'omitnan');
PQ_sink_2b   = mean(NumSink(:,idx_2B),2,'omitnan');

% Primary quadratic target uses the same boxcox accuracy convention as the
% previous PQ-behavior analysis.
Acc_2b_boxcox = mean(Acc_block_boxcox(:,idx_2B),2,'omitnan');
RT_2b = mean(RT_block(:,idx_2B),2,'omitnan');
Acc_2b_raw = mean(Acc_raw(:,5:8),2,'omitnan');
RT_2b_raw  = mean(RT(:,5:8),2,'omitnan');

stats_source_acc = fit_quadratic_pq_behavior(PQ_source_2b, Acc_2b_boxcox, sub_age, sub_gender_bi);
stats_sink_acc   = fit_quadratic_pq_behavior(PQ_sink_2b,   Acc_2b_boxcox, sub_age, sub_gender_bi);

vertex_source = stats_source_acc.vertex_pq;
vertex_sink   = stats_sink_acc.vertex_pq;

sd_source = std(PQ_source_2b,'omitnan');
sd_sink   = std(PQ_sink_2b,'omitnan');
source_signed = (PQ_source_2b - vertex_source) ./ sd_source;
sink_signed   = (PQ_sink_2b   - vertex_sink)   ./ sd_sink;
DC_signed_arm = mean([source_signed sink_signed],2,'omitnan');
DC_abs_dist   = mean([abs(source_signed) abs(sink_signed)],2,'omitnan');

arm_group = strings(numel(DC_signed_arm),1);
arm_group(DC_signed_arm < -ARM_DELTA_Z) = "left_low_PQ";
arm_group(abs(DC_signed_arm) <= ARM_DELTA_Z) = "near_optimal";
arm_group(DC_signed_arm > ARM_DELTA_Z) = "right_high_PQ";

fprintf('\n2-back D/C optimum vertices from Acc curve:\n');
fprintf('  Source/diverging vertex PQ = %.4f\n', vertex_source);
fprintf('  Sink/converging vertex PQ   = %.4f\n', vertex_sink);
fprintf('Group counts: left=%d, optimal=%d, right=%d\n', ...
    nnz(arm_group=="left_low_PQ"), nnz(arm_group=="near_optimal"), nnz(arm_group=="right_high_PQ"));

%% ---------------- B. Three-group core-trait standardized profiles ----------------
[ThreeGroupTable, ThreeGroupLong] = compute_three_group_core_traits(behav_info, core_trait_specs, ...
    arm_group, sub_age, sub_gender_bi, Acc_2b_raw, RT_2b_raw);

%% ---------------- Save all outputs ----------------
out_prefix = fullfile(cfg.outputDir,'behaviour_group_profiles');

SubjectMetricTable = table(subj_num(:), string(SubjectList(:)), PQ_source_2b(:), PQ_sink_2b(:), ...
    source_signed(:), sink_signed(:), DC_signed_arm(:), DC_abs_dist(:), arm_group(:), ...
    Acc_2b_raw(:), Acc_2b_boxcox(:), RT_2b_raw(:), RT_2b(:), ...
    'VariableNames', {'Subject','SubjectList','PQ_source_2back','PQ_sink_2back', ...
    'Source_signedFromOpt','Sink_signedFromOpt','DC_signedArm','DC_absDist','ArmGroup', ...
    'Acc_2back_raw','Acc_2back_boxcox','RT_2back_raw','RT_2back_TRexpanded'});

writetable(SubjectMetricTable, [out_prefix '_subject_metrics.csv']);
writetable(ThreeGroupTable, [out_prefix '_threegroup_coretraits.csv']);
writetable(ThreeGroupLong, [out_prefix '_threegroup_coretraits_long.csv']);

save([out_prefix '_profiles.mat'], 'SubjectMetricTable','ThreeGroupTable','ThreeGroupLong', ...
    'stats_source_acc','stats_sink_acc','vertex_source','vertex_sink','ARM_DELTA_Z');

outputs=struct('prefix',out_prefix,'subjects',numel(subj_num),'source_vertex',vertex_source,'sink_vertex',vertex_sink);
end

function core_specs = get_core_plot_specs()
% {Source, VariableName, DisplayLabel, HigherIsBetter}
% Source = 'computed' uses values computed from the fMRI WM task files.
% Source = 'csv' uses HCP unrestricted CSV columns.
core_specs = {
    'computed','Acc_2back_raw','WM 2-back accuracy',true;
    'computed','RT_2back_raw','WM 2-back speed',false;
    'csv','CogFluidComp_Unadj','Fluid cognition',true;
    'csv','CogTotalComp_Unadj','Total cognition',true;
    'csv','PMAT24_A_CR','Matrix reasoning',true;
    'csv','ListSort_Unadj','List sorting WM',true;
    'csv','ProcSpeed_Unadj','Processing speed',true;
    'csv','PercStress_Unadj','Low perceived stress',false;
    };
end

function age_mid = parse_hcp_age(age_col)
    if iscell(age_col)
        age_cell = age_col;
    else
        age_cell = cellstr(age_col);
    end
    age_mid = nan(numel(age_cell),1);
    for i = 1:numel(age_cell)
        s = char(age_cell{i});
        parts = split(s,'-');
        if numel(parts) >= 2
            age_mid(i) = (str2double(parts{1}) + str2double(parts{2})) / 2;
        else
            age_mid(i) = str2double(s);
        end
    end
end

function x = table_column_to_numeric(col)
    if isnumeric(col) || islogical(col)
        x = double(col);
    elseif iscell(col)
        x = str2double(col);
    elseif isstring(col) || ischar(col) || iscategorical(col)
        x = str2double(string(col));
    else
        try
            x = double(col);
        catch
            x = nan(size(col));
        end
    end
    x = x(:);
end

function z = zscore_safe(x)
    z = nan(size(x));
    valid = isfinite(x);
    if nnz(valid) < 3
        return;
    end
    mu = mean(x(valid),'omitnan');
    sd = std(x(valid),'omitnan');
    if ~isfinite(sd) || sd <= eps
        return;
    end
    z(valid) = (x(valid)-mu)./sd;
end

function r = residualize_covariates(y, age, gender)
    r = nan(size(y));
    valid = isfinite(y) & isfinite(age) & isfinite(gender);
    if nnz(valid) < 10
        return;
    end
    tbl = table(y(valid), age(valid), gender(valid), 'VariableNames', {'Y','Age','Gender'});
    try
        mdl = fitlm(tbl, 'Y ~ Age + Gender');
        r(valid) = mdl.Residuals.Raw;
    catch
    end
end

function stats = fit_quadratic_pq_behavior(x_pq, y_behav, age, gender)
    stats = struct('is_valid',false,'n',0,'beta_pq',nan,'beta_pq2',nan, ...
        'p_pq',nan,'p_quad',nan,'adjR2_linear',nan,'adjR2_quad',nan, ...
        'delta_adjR2',nan,'vertex_z',nan,'vertex_pq',nan, ...
        'vertex_in_range',false,'pq_min',nan,'pq_max',nan,'shape','NA');

    x_pq = x_pq(:); y_behav = y_behav(:); age = age(:); gender = gender(:);
    valid = isfinite(x_pq) & isfinite(y_behav) & isfinite(age) & isfinite(gender);
    stats.n = sum(valid);
    if stats.n < 10; return; end

    x = x_pq(valid); y = y_behav(valid); a = age(valid); g = gender(valid);
    x_mu = mean(x,'omitnan'); x_sd = std(x,'omitnan');
    if ~isfinite(x_sd) || x_sd <= eps; return; end
    xz = (x - x_mu) ./ x_sd;
    tbl = table(y, xz, xz.^2, a, g, 'VariableNames', {'Y','PQz','PQz2','Age','Gender'});
    try
        mdl_linear = fitlm(tbl, 'Y ~ PQz + Age + Gender');
        mdl_quad = fitlm(tbl, 'Y ~ PQz + PQz2 + Age + Gender');
    catch
        return;
    end
    row = mdl_quad.Coefficients.Properties.RowNames;
    b1 = mdl_quad.Coefficients.Estimate(strcmp(row,'PQz'));
    b2 = mdl_quad.Coefficients.Estimate(strcmp(row,'PQz2'));
    stats.beta_pq = b1;
    stats.beta_pq2 = b2;
    stats.p_pq = mdl_quad.Coefficients.pValue(strcmp(row,'PQz'));
    stats.p_quad = mdl_quad.Coefficients.pValue(strcmp(row,'PQz2'));
    stats.adjR2_linear = mdl_linear.Rsquared.Adjusted;
    stats.adjR2_quad = mdl_quad.Rsquared.Adjusted;
    stats.delta_adjR2 = stats.adjR2_quad - stats.adjR2_linear;
    if isfinite(b1) && isfinite(b2) && abs(b2) > eps
        stats.vertex_z = -b1 ./ (2*b2);
        stats.vertex_pq = x_mu + x_sd .* stats.vertex_z;
    end
    stats.pq_min = min(x); stats.pq_max = max(x);
    stats.vertex_in_range = isfinite(stats.vertex_pq) && stats.vertex_pq >= stats.pq_min && stats.vertex_pq <= stats.pq_max;
    if b2 < 0
        stats.shape = 'inverted-U / peak';
    elseif b2 > 0
        stats.shape = 'U-shaped / trough';
    else
        stats.shape = 'flat';
    end
    stats.is_valid = true;
end

function [ThreeGroupTable, ThreeGroupLong] = compute_three_group_core_traits(behav_info, core_specs, arm_group, age, gender, Acc_2b_raw, RT_2b_raw)
    ThreeGroupTable = table();
    ThreeGroupLong = table();


    for i = 1:size(core_specs,1)
        source = core_specs{i,1}; varname = core_specs{i,2}; label = core_specs{i,3}; higher_is_better = core_specs{i,4};
        if strcmp(source,'computed')
            if strcmp(varname,'Acc_2back_raw')
                y = Acc_2b_raw;
            elseif strcmp(varname,'RT_2back_raw')
                y = RT_2b_raw;
            else
                continue;
            end
        else
            if ~ismember(varname, behav_info.Properties.VariableNames)
                continue;
            end
            y = table_column_to_numeric(behav_info.(varname));
        end
        y_plot = y;
        if ~higher_is_better
            y_plot = -y_plot;
        end
        yr = residualize_covariates(y_plot, age, gender);
        yz = zscore_safe(yr);
        L = arm_group=="left_low_PQ" & isfinite(yz);
        O = arm_group=="near_optimal" & isfinite(yz);
        R = arm_group=="right_high_PQ" & isfinite(yz);
        if nnz(L)<5 || nnz(O)<5 || nnz(R)<5
            continue;
        end
        means = [mean(yz(L),'omitnan'), mean(yz(O),'omitnan'), mean(yz(R),'omitnan')];
        sems  = [std(yz(L),'omitnan')/sqrt(nnz(L)), std(yz(O),'omitnan')/sqrt(nnz(O)), std(yz(R),'omitnan')/sqrt(nnz(R))];
        newrow = table({source}, {varname}, {label}, logical(higher_is_better), nnz(L), nnz(O), nnz(R), ...
            means(1), means(2), means(3), sems(1), sems(2), sems(3), ...
            'VariableNames', {'Source','Variable','Label','HigherIsBetter','N_left','N_opt','N_right', ...
            'MeanZ_left','MeanZ_opt','MeanZ_right','SEM_left','SEM_opt','SEM_right'});
        ThreeGroupTable = [ThreeGroupTable; newrow]; %#ok<AGROW>

        groups = {'left/low-PQ'; 'near-optimal'; 'right/high-PQ'};
        for g = 1:3
            newlong = table({label}, groups(g), means(g), sems(g), ...
                'VariableNames', {'Label','Group','MeanZ','SEM'});
            ThreeGroupLong = [ThreeGroupLong; newlong]; %#ok<AGROW>
        end
    end
end
