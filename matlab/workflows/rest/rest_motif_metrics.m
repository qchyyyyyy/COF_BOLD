function rest_motif_metrics(save_method,i_batch,batch_size,cfg)
%STAT_VVF_HCP_YA Compute participant- and group-level HCP-YA motif metrics.
%
% Usage for regenerated REST1 LR or task-state velocity fields:
%   rest_motif_metrics('participant', 1, 299, cfg)
%   rest_motif_metrics('group', [], [], cfg)
%
% The third participant argument is the number of subjects per batch when
% cfg.batchMode is 'contiguous'. When cfg.batchMode is 'balanced', it is the
% total number of batches; the subjects are split into contiguous batches
% whose sizes differ by at most one. If cfg is omitted, historical calls with
% i_batch and batch_size retain the original strided batching semantics.

if nargin < 1 || isempty(save_method)
    save_method = 'participant';
end
if isstring(save_method)
    save_method = char(save_method);
end
if ~ischar(save_method) || ~ismember(save_method, {'participant', 'group'})
    error('save_method must be ''participant'' or ''group''.');
end
if nargin < 2
    i_batch = [];
end
if nargin < 3
    batch_size = [];
end
cfgWasSupplied = nargin >= 4 && ~isempty(cfg);
if nargin < 4 || isempty(cfg)
    cfg = struct();
end
if ~isstruct(cfg)
    error('cfg must be a struct.');
end

% Preserve the historical HCPex regional labels and binary motif threshold.
% The shared configuration supplies REST1 LR defaults while allowing the
% historical task-state sessions to be selected explicitly.
if ~isfield(cfg, 'parcellation') || isempty(cfg.parcellation)
    cfg.parcellation = 'HCPex_22';
end
if ~isfield(cfg, 'do_binarize') || isempty(cfg.do_binarize)
    cfg.do_binarize = 1;
end
if ~isfield(cfg, 'n_std') || isempty(cfg.n_std)
    cfg.n_std = 2;
end
if ~isfield(cfg, 'batchMode') || isempty(cfg.batchMode)
    if ~cfgWasSupplied && ~isempty(i_batch) && ~isempty(batch_size)
        cfg.batchMode = 'strided';
    else
        cfg.batchMode = 'contiguous';
    end
end
userSetVvfDataDir = isfield(cfg, 'vvfDataDir') && ~isempty(cfg.vvfDataDir);

cfg = rest_config(cfg);
if ~userSetVvfDataDir
    if contains(cfg.Type, 'rfMRI')
        cfg.vvfDataDir = fullfile(cfg.dataRoot, 'processed', cfg.Site, cfg.Experiment);
    else
        cfg.vvfDataDir = fullfile(cfg.dataRoot, cfg.Site, cfg.Experiment);
    end
end
if ~isfield(cfg, 'saveSubjectList') || isempty(cfg.saveSubjectList)
    % Participant jobs may run concurrently. Only the group job writes the
    % resolved list by default, avoiding races and cohort clobbering.
    cfg.saveSubjectList = strcmp(save_method, 'group');
end


parcellation = cfg.parcellation;
hemisphere = 'LR';
[~,label_L,~,~,V_mask_L] = load_cortical_atlas(parcellation,'L',cfg.dataRoot);
[~,label_R,~,~,V_mask_R] = load_cortical_atlas(parcellation,'R',cfg.dataRoot);

if strcmp(parcellation,'HCPex_22')
    label_R(label_R>0) = label_R(label_R>0) + 22;
end

V_mask = [V_mask_L V_mask_R];
V_label_L = max(label_L(1:end-1,1:end-1),max(label_L(2:end,1:end-1),max(label_L(1:end-1,2:end),label_L(2:end,2:end))));
V_label_L(V_mask_L==0)=0;
V_label_R = max(label_R(1:end-1,1:end-1),max(label_R(2:end,1:end-1),max(label_R(1:end-1,2:end),label_R(2:end,2:end))));
V_label_R(V_mask_R==0)=0;
V_label = [V_label_L V_label_R];
clearvars  label_L label_R  V_mask_L V_mask_R V_label_L V_label_R
num_regions = max(V_label, [], 'all');

Methods = normalize_methods(cfg.Method);
Site = cfg.Site;
Type = cfg.Type;
PE = cfg.PE;

Session = cfg.Session;
if ischar(Session)
    Session = {Session};
elseif isstring(Session)
    Session = cellstr(Session);
end

for i_session = 1:numel(Session)


    Experiment = [Type '_' Session{i_session} '_' PE];

    if numel(Session) == 1
        data_dir = cfg.vvfDataDir;
        participant_dir = cfg.participantMetricDir;
        group_dir = cfg.groupMetricDir;
    else
        data_dir = fullfile(cfg.dataRoot, Site, Experiment);
        results_dir = fullfile(cfg.dataRoot, 'results', Site, Experiment);
        participant_dir = fullfile(results_dir, 'participant');
        group_dir = fullfile(results_dir, 'group');
    end

    if ~exist(data_dir, 'dir')
        error('Velocity-field data directory not found: %s', data_dir);
    end

    [~,~] = mkdir(participant_dir);
    [~,~] = mkdir(group_dir);

    % Prefer the EXP-005 manifest when supplied, otherwise use a saved list
    % or discover subjects with all requested left/right velocity fields.
    [subjectList, sourceSubjectList] = resolve_subject_list( ...
        cfg, data_dir, group_dir, Experiment, Methods);
    if isempty(subjectList)
        error('No subjects with all requested velocity fields found under: %s', data_dir);
    end

    num_subjects = length(subjectList);
    c = [0; 1; -1]; % 3*1

    if strcmp(save_method, 'participant')
        if ~contains(Experiment, 'rfMRI')
            [Events, ~] = load_task_events( ...
                group_dir, Experiment, subjectList, sourceSubjectList);

            EV_Fixation = Events; EV_Fixation(Events~=0) = 0; EV_Fixation(Events==0) = 1;

            if strcmp(Session{i_session}, 'MOTOR')
                % HCP MOTOR encodes fixation as 0 and task blocks as 1.
                taskEventCode = 1;
            else
                taskEventCode = 2;
                EV_Control = Events; EV_Control(Events~=1) = 0; EV_Control(Events==1) = 1;
            end

            EV_Task = Events; EV_Task(Events~=taskEventCode) = 0;
            EV_Task(Events==taskEventCode) = 1;

            if strcmp(Session{i_session}, 'MOTOR')
                Comparison = {'Task', 'Fixation'};
            else
                Comparison = {'Control', 'Fixation'; 'Task', 'Fixation'; 'Task', 'Control'};
            end
        end

        if cfg.saveSubjectList
            save(fullfile(group_dir, [Experiment '_SubList.mat']), 'subjectList');
        end

        idx_sub = select_batch_indices(i_batch, batch_size, num_subjects, cfg.batchMode);

        for i_method = length(Methods):-1:1

            for i_sub = idx_sub
                sub_id = subjectList{i_sub};
                disp(['Task: ', Session{i_session}, ' | Sub: ', num2str(i_sub) '/' num2str(num_subjects)]);

                sub_out_dir = fullfile(participant_dir, sub_id);
                [~,~] = mkdir(sub_out_dir);

                sub_out_dir = fullfile(participant_dir, sub_id);
                if ~exist(sub_out_dir, 'dir'), mkdir(sub_out_dir); end

                if participant_outputs_complete(sub_out_dir, sub_id, Experiment, hemisphere, Methods{i_method}, parcellation)
                    continue
                end

                % Load Velocity Fields
                task_dir = fullfile(data_dir, sub_id);
                L = load(fullfile(task_dir, ['VelField_lh_', Methods{i_method}, '.mat']), 'VelField');
                R = load(fullfile(task_dir, ['VelField_rh_', Methods{i_method}, '.mat']), 'VelField');

                if strcmp(Methods{i_method}, 'PG')
                    Ux = cat(2, L.VelField.Ux(:, 1:end-1, :), R.VelField.Ux);
                    Uy = cat(2, L.VelField.Uy(:, 1:end-1, :), R.VelField.Uy);
                    Ux = convn(Ux, ones(2,2,2)/8, 'valid');
                    Uy = convn(Uy, ones(2,2,2)/8, 'valid');
                else
                    Ux = cat(2, L.VelField.Ux, R.VelField.Ux);
                    Uy = cat(2, L.VelField.Uy, R.VelField.Uy);
                end
                clearvars L R;

                % Calculate Divergence and Curl
                [DivMap, CurlMap] = velocity_derivatives(Ux, Uy, V_mask, ...
                    cfg.filtSigma, cfg.do_zscore, cfg.do_binarize, cfg.n_std);

                % Calculate Order Parameter
                T = size(DivMap, 3);
                if ~contains(Experiment, 'rfMRI') && size(Events, 2) < T
                    error(['Task event matrix for %s has %d time points, but ' ...
                        'the velocity field has %d frames.'], Experiment, ...
                        size(Events, 2), T);
                end
                U_mag = sqrt(Ux.^2 + Uy.^2);
                Ux = Ux ./ U_mag;
                Uy = Uy ./ U_mag;
                Ux(repmat(V_mask, 1, 1, T) == 0) = nan;
                Uy(repmat(V_mask, 1, 1, T) == 0) = nan;

                imageFilter = fspecial('average', 9);
                for t = 1:T
                    Ux(:,:,t) = nanconv(Ux(:,:,t), imageFilter);
                    Uy(:,:,t) = nanconv(Uy(:,:,t), imageFilter);
                end

                order_param = sqrt(Ux.^2 + Uy.^2);
                order_param(repmat(V_mask, 1, 1, T) == 0) = nan;
                clearvars Ux Uy U_mag;

                % --- 1 Save Order Parameter ---
                % Global
                order_param_global = squeeze(mean(order_param, [1 2], 'omitnan'));
                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_global', '-v7.3');

                % Regional
                order_param_region = zeros(num_regions, T);
                for i_region = 1:num_regions
                    mask_rep = repmat(V_label == i_region, 1, 1, T);
                    tmp_param = order_param;
                    tmp_param(~mask_rep) = nan;
                    order_param_region(i_region, :) = squeeze(mean(tmp_param, [1 2], 'omitnan'));
                end
                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_region', '-v7.3');

                % Voxel-wise
                if contains(Experiment, 'rfMRI')
                    order_param_mean = mean(order_param, 3, 'omitnan');
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_mean', '-v7.3');
                else
                    if strcmp(Session{i_session},'MOTOR')
                        order_param_Fixation = mean(order_param(:, :, Events(i_sub, 1:T) == 0), 3, 'omitnan');
                        order_param_Task = mean(order_param(:, :, Events(i_sub, 1:T) == taskEventCode), 3, 'omitnan');
                        save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_Fixation', 'order_param_Task', '-v7.3');
                    else
                        order_param_Fixation = mean(order_param(:, :, Events(i_sub, 1:T) == 0), 3, 'omitnan');
                        order_param_Control = mean(order_param(:, :, Events(i_sub, 1:T) == 1), 3, 'omitnan');
                        order_param_Task = mean(order_param(:, :, Events(i_sub, 1:T) == taskEventCode), 3, 'omitnan');
                        save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_Fixation', 'order_param_Control', 'order_param_Task', '-v7.3');
                    end
                    
                    % effect size
                    data = grid_to_vector(order_param,V_mask)'; % T*V
                    cope = []; 
                    for i_compar = size(Comparison,1):-1:1
                        eval(['design_mat = [ones(size(Events(i_sub, 1:T)'')) EV_' Comparison{i_compar,1} '(i_sub,1:T)'' EV_' Comparison{i_compar,2} '(i_sub,1:T)''];']);
                        cope(:,:,i_compar) = fit_gls_ar1(data, design_mat, c);
                    end
                    save(fullfile(sub_out_dir, [sub_id '_' Experiment '_orderParam_level1_' hemisphere '_' Methods{i_method} '.mat']),...
                        'cope');
                end

                % --- 2 Save Pattern Sizes ---
                SizeSource = zeros(1, T); SizeSink = zeros(1, T);
                SizeCCVortex = zeros(1, T); SizeCVortex = zeros(1, T);

                SizeSource_region = zeros(num_regions, T); SizeSink_region = zeros(num_regions, T);
                SizeCCVortex_region = zeros(num_regions, T); SizeCVortex_region = zeros(num_regions, T);

                for t = 1:T
                    % Global size
                    [SizeSource(1,t), SizeSink(1,t)] = measure_motif_sizes(DivMap(:,:,t));
                    [SizeCCVortex(1,t), SizeCVortex(1,t)] = measure_motif_sizes(CurlMap(:,:,t));

                    % Regional size
                    for i_region = 1:num_regions
                        region_mask = (V_label == i_region);

                        tmp_DivMap = DivMap(:,:,t); tmp_DivMap(~region_mask) = 0;
                        tmp_CurlMap = CurlMap(:,:,t); tmp_CurlMap(~region_mask) = 0;

                        [SizeSource_region(i_region, t), SizeSink_region(i_region, t)] = measure_motif_sizes(tmp_DivMap);
                        [SizeCCVortex_region(i_region, t), SizeCVortex_region(i_region, t)] = measure_motif_sizes(tmp_CurlMap);
                    end
                end

                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                    'SizeSource', 'SizeSink', 'SizeCCVortex', 'SizeCVortex', '-v7.3');

                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                    'SizeSource_region', 'SizeSink_region', 'SizeCCVortex_region', 'SizeCVortex_region', '-v7.3');

                % --- 3 Save Pattern Quantity ---
                NumSource = zeros(1, T); NumSink = zeros(1, T);
                NumCCVortex = zeros(1, T); NumCVortex = zeros(1, T);

                NumSource_region = zeros(num_regions, T); NumSink_region = zeros(num_regions, T);
                NumCCVortex_region = zeros(num_regions, T); NumCVortex_region = zeros(num_regions, T);

                for t = 1:T
                    % Global size
                    [NumSource(1,t), NumSink(1,t)] = count_motifs(DivMap(:,:,t));
                    [NumCCVortex(1,t), NumCVortex(1,t)] = count_motifs(CurlMap(:,:,t));

                    % Regional size
                    for i_region = 1:num_regions
                        region_mask = (V_label == i_region);

                        tmp_DivMap = DivMap(:,:,t); tmp_DivMap(~region_mask) = 0;
                        tmp_CurlMap = CurlMap(:,:,t); tmp_CurlMap(~region_mask) = 0;

                        [NumSource_region(i_region, t), NumSink_region(i_region, t)] = count_motifs(tmp_DivMap);
                        [NumCCVortex_region(i_region, t), NumCVortex_region(i_region, t)] = count_motifs(tmp_CurlMap);
                    end
                end

                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                    'NumSource', 'NumSink', 'NumCCVortex', 'NumCVortex', '-v7.3');

                save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                    'NumSource_region', 'NumSink_region', 'NumCCVortex_region', 'NumCVortex_region', '-v7.3');

                % --- 3 Save voxel-wise Pattern density ---


                SoM = DivMap;SoM(SoM<0)=0;
                SiM = -DivMap;SiM(SiM<0)=0;
                CCM = CurlMap;CCM(CCM<0)=0;
                CM = -CurlMap;CM(CM<0)=0;

                if contains(Experiment,'rfMRI')
                    SourceMap = mean(SoM,3,'omitnan');
                    SinkMap = mean(SiM,3,'omitnan');
                    CCVortexMap = mean(CCM,3,'omitnan');
                    CVortexMap = mean(CM,3,'omitnan');

                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                        'SourceMap','SinkMap','CCVortexMap','CVortexMap','-v7.3');
                else
                    SourceMap_Fixation = mean(SoM(:,:,Events(i_sub,1:T)==0),3,'omitnan');
                    SinkMap_Fixation = mean(SiM(:,:,Events(i_sub,1:T)==0),3,'omitnan');
                    CCVortexMap_Fixation = mean(CCM(:,:,Events(i_sub,1:T)==0),3,'omitnan');
                    CVortexMap_Fixation = mean(CM(:,:,Events(i_sub,1:T)==0),3,'omitnan');

                    if ~strcmp(Session{i_session},'MOTOR')
                        SourceMap_Control = mean(SoM(:,:,Events(i_sub,1:T)==1),3,'omitnan');
                        SinkMap_Control = mean(SiM(:,:,Events(i_sub,1:T)==1),3,'omitnan');
                        CCVortexMap_Control = mean(CCM(:,:,Events(i_sub,1:T)==1),3,'omitnan');
                        CVortexMap_Control = mean(CM(:,:,Events(i_sub,1:T)==1),3,'omitnan');
                    end

                    SourceMap_Task = mean(SoM(:,:,Events(i_sub,1:T)==taskEventCode),3,'omitnan');
                    SinkMap_Task = mean(SiM(:,:,Events(i_sub,1:T)==taskEventCode),3,'omitnan');
                    CCVortexMap_Task = mean(CCM(:,:,Events(i_sub,1:T)==taskEventCode),3,'omitnan');
                    CVortexMap_Task = mean(CM(:,:,Events(i_sub,1:T)==taskEventCode),3,'omitnan');

                    if strcmp(Session{i_session},'MOTOR')
                        save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                            'SourceMap_Fixation', 'SinkMap_Fixation', 'CCVortexMap_Fixation', 'CVortexMap_Fixation',...
                            'SourceMap_Task','SinkMap_Task','CCVortexMap_Task','CVortexMap_Task','-v7.3');
                    else
                        save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                            'SourceMap_Fixation', 'SinkMap_Fixation', 'CCVortexMap_Fixation', 'CVortexMap_Fixation',...
                            'SourceMap_Control','SinkMap_Control','CCVortexMap_Control','CVortexMap_Control',...
                            'SourceMap_Task','SinkMap_Task','CCVortexMap_Task','CVortexMap_Task','-v7.3');
                    end

                    % effect size
                    cope = [];
                    
                    Map = {'SoM','SiM','CCM','CM'};
                    for i_pattern = 1:numel(Map)
                        eval(['data = grid_to_vector(' Map{i_pattern} ',V_mask)'';']);
                        for i_compar = size(Comparison,1):-1:1
                            eval(['design_mat = [ones(size(Events(i_sub,1:T)'')) EV_' Comparison{i_compar,1} '(i_sub,1:T)'' EV_' Comparison{i_compar,2} '(i_sub,1:T)''];']);
                            cope(:,:,i_compar,i_pattern) = fit_gls_ar1(data, design_mat, c);
                        end
                    end
                    save(fullfile(sub_out_dir, [sub_id '_' Experiment '_numPattern_level1_' hemisphere '_' Methods{i_method} '.mat']),...
                        'cope');
                end
            end
        end

    elseif strcmp(save_method, 'group')

        if cfg.saveSubjectList
            save(fullfile(group_dir, [Experiment '_SubList.mat']), 'subjectList');
        end

        % Group aggregation must never silently insert zeros for an absent
        % participant output. The participant stage is therefore a hard
        % prerequisite, including all four motif-density maps used by Fig. 1e.
        assert_participant_outputs_complete(subjectList, participant_dir, ...
            Experiment, hemisphere, Methods, parcellation);
        
        if strcmp(Session{i_session}, 'MOTOR')
            Comparison = {'Task', 'Fixation'};
        else
            Comparison = {'Control', 'Fixation'; 'Task', 'Fixation'; 'Task', 'Control'};
        end

        for i_method = length(Methods):-1:1

            for i_sub = 1:num_subjects

                sub_id = subjectList{i_sub};

                disp(['Task: ', Session{i_session}, ' | Sub: ', num2str(i_sub) '/' num2str(num_subjects)]);

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);

                    if contains(Experiment,'rfMRI')
                        if i_sub==1
                            order_param = zeros([size(V_mask) num_subjects]);
                        end
                        order_param(:,:,i_sub) = tmp.order_param_mean;
                    else
                        if i_sub==1
                            if ~strcmp(Session{i_session},'MOTOR')
                                order_param_Control = zeros([size(V_mask) num_subjects]);
                            end
                            order_param_Task = zeros([size(V_mask) num_subjects]);
                            order_param_Fixation = zeros([size(V_mask) num_subjects]);
                        end
                        if ~strcmp(Session{i_session},'MOTOR')
                            order_param_Control(:,:,i_sub) = tmp.order_param_Control;
                        end
                        order_param_Task(:,:,i_sub) = tmp.order_param_Task;
                        order_param_Fixation(:,:,i_sub) = tmp.order_param_Fixation;
                    end
                else
                    if ~contains(Experiment, 'rfMRI')
                        warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                    end
                end
                
                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_level1_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);

                    if ~contains(Experiment,'rfMRI')
                        if i_sub==1
                            order_param_cope = zeros([num_subjects nnz(V_mask==1) numel(Comparison)]);
                        end
                        order_param_cope(i_sub,:,:) = tmp.cope;
                    end
                else
                    if ~contains(Experiment, 'rfMRI')
                        warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                    end
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = length(tmp.order_param_global);
                    if i_sub==1
                        order_param_global = zeros(num_subjects, T);
                    end
                    order_param_global(i_sub,:) = tmp.order_param_global;
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = size(tmp.order_param_region, 2);
                    if i_sub==1
                        order_param_region   = zeros(num_regions, T, num_subjects);
                    end

                    order_param_region(:,:,i_sub)   = tmp.order_param_region;

                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end



                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = size(tmp.NumSource, 2);
                    if i_sub==1
                        NumSource   = zeros(num_subjects, T);
                        NumSink     = zeros(num_subjects, T);
                        NumCCVortex = zeros(num_subjects, T);
                        NumCVortex  = zeros(num_subjects, T);
                    end

                    NumSource(i_sub,:)   = tmp.NumSource;
                    NumSink(i_sub,:)     = tmp.NumSink;
                    NumCCVortex(i_sub,:) = tmp.NumCCVortex;
                    NumCVortex(i_sub,:)  = tmp.NumCVortex;
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = size(tmp.NumSource_region, 2);
                    if i_sub==1
                        NumSource_region   = zeros(num_regions, T, num_subjects);
                        NumSink_region     = zeros(num_regions, T, num_subjects);
                        NumCCVortex_region = zeros(num_regions, T, num_subjects);
                        NumCVortex_region  = zeros(num_regions, T, num_subjects);
                    end

                    NumSource_region(:,:,i_sub)   = tmp.NumSource_region;
                    NumSink_region(:,:,i_sub)     = tmp.NumSink_region;
                    NumCCVortex_region(:,:,i_sub) = tmp.NumCCVortex_region;
                    NumCVortex_region(:,:,i_sub)  = tmp.NumCVortex_region;
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = size(tmp.SizeSource, 2);
                    if i_sub==1
                        SizeSource   = zeros(num_subjects, T);
                        SizeSink     = zeros(num_subjects, T);
                        SizeCCVortex = zeros(num_subjects, T);
                        SizeCVortex  = zeros(num_subjects, T);
                    end

                    SizeSource(i_sub,:)   = tmp.SizeSource;
                    SizeSink(i_sub,:)     = tmp.SizeSink;
                    SizeCCVortex(i_sub,:) = tmp.SizeCCVortex;
                    SizeCVortex(i_sub,:)  = tmp.SizeCVortex;
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);
                    T = size(tmp.SizeSource_region, 2);
                    if i_sub==1
                        SizeSource_region   = zeros(num_regions, T, num_subjects);
                        SizeSink_region     = zeros(num_regions, T, num_subjects);
                        SizeCCVortex_region = zeros(num_regions, T, num_subjects);
                        SizeCVortex_region  = zeros(num_regions, T, num_subjects);
                    end

                    SizeSource_region(:,:,i_sub)   = tmp.SizeSource_region;
                    SizeSink_region(:,:,i_sub)     = tmp.SizeSink_region;
                    SizeCCVortex_region(:,:,i_sub) = tmp.SizeCCVortex_region;
                    SizeCVortex_region(:,:,i_sub)  = tmp.SizeCVortex_region;
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);

                    if contains(Experiment,'rfMRI')
                        if i_sub==1
                            SourceMap = zeros([size(V_mask) num_subjects]);
                            SinkMap = zeros([size(V_mask) num_subjects]);
                            CCVortexMap = zeros([size(V_mask) num_subjects]);
                            CVortexMap = zeros([size(V_mask) num_subjects]);
                        end
                        SourceMap(:,:,i_sub) = tmp.SourceMap;
                        SinkMap(:,:,i_sub) = tmp.SinkMap;
                        CCVortexMap(:,:,i_sub) = tmp.CCVortexMap;
                        CVortexMap(:,:,i_sub) = tmp.CVortexMap;

                    else
                        if i_sub==1
                            SourceMap_Fixation = zeros([size(V_mask) num_subjects]);
                            SinkMap_Fixation = zeros([size(V_mask) num_subjects]);
                            CCVortexMap_Fixation = zeros([size(V_mask) num_subjects]);
                            CVortexMap_Fixation = zeros([size(V_mask) num_subjects]);

                            if ~strcmp(Session{i_session},'MOTOR')
                                SourceMap_Control = zeros([size(V_mask) num_subjects]);
                                SinkMap_Control = zeros([size(V_mask) num_subjects]);
                                CCVortexMap_Control = zeros([size(V_mask) num_subjects]);
                                CVortexMap_Control = zeros([size(V_mask) num_subjects]);
                            end

                            SourceMap_Task = zeros([size(V_mask) num_subjects]);
                            SinkMap_Task = zeros([size(V_mask) num_subjects]);
                            CCVortexMap_Task = zeros([size(V_mask) num_subjects]);
                            CVortexMap_Task = zeros([size(V_mask) num_subjects]);

                        end
                        SourceMap_Fixation(:,:,i_sub) = tmp.SourceMap_Fixation;
                        SinkMap_Fixation(:,:,i_sub) = tmp.SinkMap_Fixation;
                        CCVortexMap_Fixation(:,:,i_sub) = tmp.CCVortexMap_Fixation;
                        CVortexMap_Fixation(:,:,i_sub) = tmp.CVortexMap_Fixation;

                        if ~strcmp(Session{i_session},'MOTOR')
                            SourceMap_Control(:,:,i_sub) = tmp.SourceMap_Control;
                            SinkMap_Control(:,:,i_sub) = tmp.SinkMap_Control;
                            CCVortexMap_Control(:,:,i_sub) = tmp.CCVortexMap_Control;
                            CVortexMap_Control(:,:,i_sub) = tmp.CVortexMap_Control;
                        end

                        SourceMap_Task(:,:,i_sub) = tmp.SourceMap_Task;
                        SinkMap_Task(:,:,i_sub) = tmp.SinkMap_Task;
                        CCVortexMap_Task(:,:,i_sub) = tmp.CCVortexMap_Task;
                        CVortexMap_Task(:,:,i_sub) = tmp.CVortexMap_Task;
                    end
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_numPattern_level1_', hemisphere, '_', Methods{i_method}, '.mat']);
                if exist(tmp_file, 'file')
                    tmp = load(tmp_file);

                    if ~contains(Experiment,'rfMRI')
                        if i_sub==1
                            pattern_cope = zeros([num_subjects nnz(V_mask==1) numel(Comparison) 4]);
                        end
                        pattern_cope(i_sub,:,:,:) = tmp.cope;
                    end
                else
                    if ~contains(Experiment, 'rfMRI')
                        warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                    end
                end

            end

            % 保存该 Method 和 Task 下的组级矩阵
            group_save_name = fullfile(group_dir, ['task-', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'order_param_global', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'order_param_region', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
            if contains(Experiment,'rfMRI')
                save(group_save_name, 'order_param', '-v7.3');
            else
                if strcmp(Session{i_session},'MOTOR')
                    save(group_save_name, 'order_param_Fixation', 'order_param_Task', '-v7.3');
                else
                    save(group_save_name, 'order_param_Fixation', 'order_param_Control', 'order_param_Task', '-v7.3');
                end
            end
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_orderParam_level1_', hemisphere, '_', Methods{i_method}, '.mat']);
            if ~contains(Experiment,'rfMRI')
                save(group_save_name, 'order_param_cope', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'NumSource', 'NumSink', 'NumCCVortex', 'NumCVortex', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'NumSource_region', 'NumSink_region', 'NumCCVortex_region', 'NumCVortex_region', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'SizeSource', 'SizeSink', 'SizeCCVortex', 'SizeCVortex', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
            save(group_save_name, 'SizeSource_region', 'SizeSink_region', 'SizeCCVortex_region', 'SizeCVortex_region', '-v7.3');
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
            if contains(Experiment,'rfMRI')
                save(group_save_name, 'SourceMap', 'SinkMap', 'CCVortexMap', 'CVortexMap', '-v7.3');
            else
                if strcmp(Session{i_session},'MOTOR')
                    save(group_save_name,...
                        'SourceMap_Fixation', 'SinkMap_Fixation', 'CCVortexMap_Fixation', 'CVortexMap_Fixation',...
                        'SourceMap_Task', 'SinkMap_Task', 'CCVortexMap_Task', 'CVortexMap_Task','-v7.3');
                else
                    save(group_save_name,...
                        'SourceMap_Fixation', 'SinkMap_Fixation', 'CCVortexMap_Fixation', 'CVortexMap_Fixation',...
                        'SourceMap_Control', 'SinkMap_Control', 'CCVortexMap_Control', 'CVortexMap_Control',...
                        'SourceMap_Task', 'SinkMap_Task', 'CCVortexMap_Task', 'CVortexMap_Task','-v7.3');
                end
            end
            disp(['Saved: ', group_save_name]);

            group_save_name = fullfile(group_dir, ['task-', Experiment, '_numPattern_level1_', hemisphere, '_', Methods{i_method}, '.mat']);
            if ~contains(Experiment,'rfMRI')
                save(group_save_name, 'pattern_cope', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end


            % (你可以用同样的逻辑，继续汇总 order_param 或 voxel-level 的数据)
        end
        disp('Group level aggregation completed!');
    end

end

end

function Methods = normalize_methods(methods)

if ischar(methods)
    Methods = {methods};
elseif isstring(methods)
    Methods = cellstr(methods(:));
elseif iscell(methods)
    Methods = methods(:);
else
    error('cfg.Method must be a character vector, string array, or cell array.');
end

for i_method = 1:numel(Methods)
    if isstring(Methods{i_method})
        Methods{i_method} = char(Methods{i_method});
    end
    if ~ischar(Methods{i_method}) || isempty(strtrim(Methods{i_method}))
        error('Each requested method must be a non-empty character vector.');
    end
    Methods{i_method} = strtrim(Methods{i_method});
end

Methods = unique(Methods, 'stable');

end

function idx_sub = select_batch_indices(i_batch, batch_size, num_subjects, batch_mode)

if isempty(i_batch)
    idx_sub = 1:num_subjects;
    return
end
if isempty(batch_size)
    % The historical two-argument form ignored an isolated i_batch and ran
    % the full cohort; retain that behavior for compatibility.
    idx_sub = 1:num_subjects;
    return
end

if ~isscalar(i_batch) || ~isfinite(i_batch) || i_batch < 1 || i_batch ~= floor(i_batch)
    error('i_batch must be a positive integer.');
end
if isempty(batch_size) || ~isscalar(batch_size) || ~isfinite(batch_size) || ...
        batch_size < 1 || batch_size ~= floor(batch_size)
    error('batch_size must be a positive integer when i_batch is supplied.');
end

if isstring(batch_mode)
    batch_mode = char(batch_mode);
end
batch_mode = lower(strtrim(batch_mode));

switch batch_mode
    case 'contiguous'
        first_idx = (i_batch - 1) * batch_size + 1;
        last_idx = min(i_batch * batch_size, num_subjects);
        if first_idx > num_subjects
            idx_sub = [];
        else
            idx_sub = first_idx:last_idx;
        end
    case 'strided'
        % Compatibility mode for the historical implementation, where
        % i_batch is the array-task index and batch_size is the number of
        % array tasks.
        idx_sub = i_batch:batch_size:num_subjects;
    case 'balanced'
        % Use the third argument as the total number of array tasks and
        % distribute the subject list as evenly as possible. For the
        % EXP-006 cohort (1,074 subjects and 299 tasks), each task therefore
        % receives either three or four subjects.
        if i_batch > batch_size
            error('i_batch %d exceeds the balanced batch count %d.', ...
                i_batch, batch_size);
        end
        first_idx = floor((i_batch - 1) * num_subjects / batch_size) + 1;
        last_idx = floor(i_batch * num_subjects / batch_size);
        if first_idx > num_subjects || last_idx < first_idx
            idx_sub = [];
        else
            idx_sub = first_idx:last_idx;
        end
    otherwise
        error('Unknown cfg.batchMode "%s". Use "contiguous", "strided", or "balanced".', batch_mode);
end

end

function [subjectList, sourceSubjectList] = resolve_subject_list( ...
        cfg, data_dir, group_dir, Experiment, Methods)

manifestFile = '';
if isfield(cfg, 'subjectManifestFile') && ~isempty(cfg.subjectManifestFile)
    manifestFile = cfg.subjectManifestFile;
end

if ~isempty(manifestFile)
    if ~exist(manifestFile, 'file')
        error('Configured subject manifest not found: %s', manifestFile);
    end
    candidateList = read_subject_manifest(manifestFile);
elseif exist(fullfile(group_dir, [Experiment '_SubList.mat']), 'file')
    saved = load(fullfile(group_dir, [Experiment '_SubList.mat']), 'subjectList');
    if ~isfield(saved, 'subjectList')
        error('Saved subject list is missing the variable subjectList: %s', ...
            fullfile(group_dir, [Experiment '_SubList.mat']));
    end
    candidateList = normalize_subject_list(saved.subjectList);
else
    candidateList = discover_subjects_with_methods(data_dir, Methods);
end

candidateList = normalize_subject_list(candidateList);
candidateList = unique(candidateList, 'stable');
sourceSubjectList = candidateList;
hasVelocityFields = false(size(candidateList));
for i_sub = 1:numel(candidateList)
    subjectDir = fullfile(data_dir, candidateList{i_sub});
    hasVelocityFields(i_sub) = exist(subjectDir, 'dir') == 7;
    for i_method = 1:numel(Methods)
        hasVelocityFields(i_sub) = hasVelocityFields(i_sub) && ...
            exist(fullfile(subjectDir, ['VelField_lh_' Methods{i_method} '.mat']), 'file') == 2 && ...
            exist(fullfile(subjectDir, ['VelField_rh_' Methods{i_method} '.mat']), 'file') == 2;
    end
end

if any(~hasVelocityFields)
    warning('rest_motif_metrics:IncompleteVelocityFields', ...
        'Excluding %d subjects without all requested left/right velocity fields.', ...
        nnz(~hasVelocityFields));
end
subjectList = candidateList(hasVelocityFields);

end

function [Events, eventSubjectList] = load_task_events( ...
        group_dir, Experiment, subjectList, sourceSubjectList)

eventsFile = fullfile(group_dir, [Experiment '_EVs.mat']);
if ~exist(eventsFile, 'file')
    error('Task event file not found: %s', eventsFile);
end

saved = load(eventsFile, 'Events');
if ~isfield(saved, 'Events')
    error('Task event file is missing the variable Events: %s', eventsFile);
end
Events = saved.Events;
if ~(isnumeric(Events) || islogical(Events)) || ~ismatrix(Events)
    error('Events must be a numeric two-dimensional matrix: %s', eventsFile);
end

% Prefer subject IDs stored with the event matrix. The historical task files
% often store only Events, so fall back to the contemporaneous SubList file
% and finally to the source list used to discover/filter velocity fields.
eventSubjectList = {};
eventInfo = whos('-file', eventsFile);
eventNames = {eventInfo.name};
if any(strcmp(eventNames, 'subjectList'))
    savedIds = load(eventsFile, 'subjectList');
    eventSubjectList = normalize_subject_list(savedIds.subjectList);
elseif exist(fullfile(group_dir, [Experiment '_SubList.mat']), 'file')
    savedIds = load(fullfile(group_dir, [Experiment '_SubList.mat']), 'subjectList');
    if isfield(savedIds, 'subjectList')
        eventSubjectList = normalize_subject_list(savedIds.subjectList);
    end
end
if isempty(eventSubjectList)
    eventSubjectList = sourceSubjectList;
end
eventSubjectList = unique(eventSubjectList, 'stable');

if numel(eventSubjectList) ~= size(Events, 1)
    error(['Task event rows (%d) do not match the associated subject list ' ...
        '(%d) for %s.'], size(Events, 1), numel(eventSubjectList), Experiment);
end

[isPresent, eventRows] = ismember(subjectList, eventSubjectList);
if any(~isPresent)
    missing = subjectList(~isPresent);
    nShown = min(numel(missing), 10);
    error(['Task event file has no row for %d velocity-field subject(s) in %s. ' ...
        'First missing IDs: %s'], numel(missing), Experiment, ...
        strjoin(missing(1:nShown), ', '));
end

% Reorder/subset event rows after velocity-field completeness filtering so
% Events(i_sub,:) always describes subjectList{i_sub}.
Events = Events(eventRows, :);

end

function subjectList = read_subject_manifest(manifestFile)

fid = fopen(manifestFile, 'r');
if fid < 0
    error('Cannot open subject manifest: %s', manifestFile);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
raw = textscan(fid, '%s');
subjectList = raw{1};
subjectList = subjectList(~cellfun(@isempty, subjectList));

end

function subjectList = discover_subjects_with_methods(data_dir, Methods)

items = dir(data_dir);
subjectList = {};
for i_item = 1:numel(items)
    if ~items(i_item).isdir || strcmp(items(i_item).name, '.') || strcmp(items(i_item).name, '..')
        continue
    end
    subjectDir = fullfile(data_dir, items(i_item).name);
    hasAllMethods = true;
    for i_method = 1:numel(Methods)
        hasAllMethods = hasAllMethods && ...
            exist(fullfile(subjectDir, ['VelField_lh_' Methods{i_method} '.mat']), 'file') == 2 && ...
            exist(fullfile(subjectDir, ['VelField_rh_' Methods{i_method} '.mat']), 'file') == 2;
    end
    if hasAllMethods
        subjectList{end+1, 1} = items(i_item).name; %#ok<AGROW>
    end
end
subjectList = sort(subjectList);

end

function subjectList = normalize_subject_list(subjectList)

if ischar(subjectList)
    subjectList = cellstr(subjectList);
elseif isstring(subjectList)
    subjectList = cellstr(subjectList(:));
elseif isnumeric(subjectList)
    subjectList = cellstr(num2str(subjectList(:)));
elseif iscell(subjectList)
    subjectList = subjectList(:);
else
    error('Subject list must be text, numeric IDs, or a cell array.');
end

for i_sub = 1:numel(subjectList)
    if isnumeric(subjectList{i_sub})
        subjectList{i_sub} = num2str(subjectList{i_sub});
    elseif isstring(subjectList{i_sub})
        subjectList{i_sub} = char(subjectList{i_sub});
    end
    if ~ischar(subjectList{i_sub})
        error('Subject IDs must be character vectors.');
    end
    subjectList{i_sub} = strtrim(subjectList{i_sub});
end
subjectList = subjectList(~cellfun(@isempty, subjectList));

end

function tf = participant_outputs_complete(sub_out_dir, sub_id, Experiment, hemisphere, method, parcellation)

requiredFiles = participant_output_files(sub_out_dir, sub_id, Experiment, hemisphere, method, parcellation);
tf = all(cellfun(@(fileName) exist(fileName, 'file') == 2, requiredFiles));

end

function assert_participant_outputs_complete(subjectList, participant_dir, ...
        Experiment, hemisphere, Methods, parcellation)

missing = {};
for i_method = 1:numel(Methods)
    for i_sub = 1:numel(subjectList)
        sub_id = subjectList{i_sub};
        requiredFiles = participant_output_files( ...
            fullfile(participant_dir, sub_id), sub_id, Experiment, hemisphere, ...
            Methods{i_method}, parcellation);
        for i_file = 1:numel(requiredFiles)
            if exist(requiredFiles{i_file}, 'file') ~= 2
                missing{end+1, 1} = requiredFiles{i_file}; %#ok<AGROW>
            end
        end
    end
end

if ~isempty(missing)
    nShown = min(numel(missing), 10);
    detail = strjoin(missing(1:nShown), sprintf('\n'));
    error('Participant outputs are incomplete: %d file(s) missing. First missing files:\n%s', ...
        numel(missing), detail);
end

end

function requiredFiles = participant_output_files(sub_out_dir, sub_id, Experiment, hemisphere, method, parcellation)

prefix = [sub_id '_' Experiment '_'];
requiredFiles = {
    fullfile(sub_out_dir, [prefix 'orderParam_global_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'orderParam_region-' parcellation '_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'orderParam_voxel_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'sizePattern_global_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'sizePattern_region-' parcellation '_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'numPattern_global_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'numPattern_region-' parcellation '_' hemisphere '_' method '.mat'])
    fullfile(sub_out_dir, [prefix 'numPattern_voxel_' hemisphere '_' method '.mat'])
    };

% REST has no event-wise first-level contrast files. Task runs retain the
% historical level-1 outputs and therefore require them for aggregation.
if ~contains(Experiment, 'rfMRI')
    requiredFiles{end+1} = fullfile(sub_out_dir, ...
        [prefix 'orderParam_level1_' hemisphere '_' method '.mat']);
    requiredFiles{end+1} = fullfile(sub_out_dir, ...
        [prefix 'numPattern_level1_' hemisphere '_' method '.mat']);
end

end
