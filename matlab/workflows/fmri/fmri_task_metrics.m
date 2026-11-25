function fmri_task_metrics(save_method,i_batch,batch_size,analysisFlags)
% Historical participant feature and COPE producer for current Fig2g/ED8.
if nargin<4, analysisFlags=struct(); end

assert(~isempty(getenv('COF_DATA_ROOT')) && ~isempty(getenv('COF_RESULTS_ROOT')) && ~isempty(getenv('COF_TASK_OUTPUT_ROOT')),'Set COF_DATA_ROOT, COF_RESULTS_ROOT and a fresh COF_TASK_OUTPUT_ROOT.');
assert(~strcmp(getenv('COF_RESULTS_ROOT'),getenv('COF_TASK_OUTPUT_ROOT')),'Input and output roots must differ.');

parcellation = 'HCPex_22';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';
[~,label_L,~,~,V_mask_L] = load_cortical_atlas(parcellation,'L');
[~,label_R,~,~,V_mask_R] = load_cortical_atlas(parcellation,'R');

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

% Methods = {'OF','PG','COF'};
Methods = {'COF'};

Site = 'HCP-YA';
Type = 'tfMRI';
PE = 'LR';

pattern_names_old = {'Source','Sink','CCVortex','CVortex'};

Session = {'EMOTION', 'GAMBLING', 'LANGUAGE', 'RELATIONAL', 'SOCIAL', 'WM'}; % , 'MOTOR'

for i_session = 1:numel(Session)

    Experiment = [Type '_' Session{i_session} '_' PE];

    data_dir = fullfile(getenv('COF_DATA_ROOT'), Site, Experiment);
    results_dir = fullfile(getenv('COF_TASK_OUTPUT_ROOT'), Site, Experiment);
    metadata_dir=fullfile(getenv('COF_RESULTS_ROOT'),Site,Experiment,'group');
    participant_dir = fullfile(results_dir, 'participant');
    group_dir = fullfile(results_dir, 'group');

    [~,~] = mkdir(participant_dir);
    [~,~] = mkdir(group_dir);

    % Access the list of subject files
    load(fullfile(metadata_dir, [Experiment '_SubList.mat']),'subjectList');

    num_subjects = length(subjectList);
    c = [0; 1; -1]; % 3*1

    cfg = [];
    cfg.do_pattern_quantity_global = false;
    cfg.do_pattern_quantity_region = false;
    cfg.do_pattern_size_global = false;
    cfg.do_pattern_size_region = false;
    cfg.do_pattern_voxel = true;
    cfg.do_order_param_global = false;
    cfg.do_order_param_region = false;
    cfg.do_order_param_voxel = false;
    names=fieldnames(analysisFlags);
    for q=1:numel(names)
        assert(isfield(cfg,names{q}),'Unknown feature flag: %s',names{q});
        cfg.(names{q})=analysisFlags.(names{q});
    end

    if strcmp(save_method, 'participant')
        load(fullfile(metadata_dir, [Experiment '_EVs.mat']), 'Events');

        event_regressors.Fixation = Events; event_regressors.Fixation(Events~=0) = 0; event_regressors.Fixation(Events==0) = 1;

        if ~strcmp(Session{i_session}, 'MOTOR')
            event_regressors.Control = Events; event_regressors.Control(Events~=1) = 0; event_regressors.Control(Events==1) = 1;
            event_regressors.Task = Events; event_regressors.Task(Events~=2) = 0; event_regressors.Task(Events==2) = 1;
        else
            n_conditions = max(Events,[],'all');
            for i_condition = 1:n_conditions
                event_regressors.(['Task' num2str(i_condition)]) = double(Events==i_condition);
            end
        end

        if strcmp(Session{i_session}, 'MOTOR')
            Comparison = {'Task1', 'Fixation'; 'Task2', 'Fixation'; 'Task3', 'Fixation'; 'Task4', 'Fixation'; 'Task5', 'Fixation'};
        else
            Comparison = {'Control', 'Fixation'; 'Task', 'Fixation'; 'Task', 'Control'};
        end

        if nargin >= 3 && ~isempty(i_batch) && ~isempty(batch_size)
            idx_sub = i_batch:batch_size:num_subjects;
        else
            idx_sub = 1:num_subjects;
        end

        for i_method = length(Methods):-1:1

            for i_sub = idx_sub
                sub_id = subjectList{i_sub};
                disp(['Task: ', Session{i_session}, ' | Sub: ', num2str(i_sub) '/' num2str(num_subjects)]);

                sub_out_dir = fullfile(participant_dir, sub_id);
                [~,~] = mkdir(sub_out_dir);

                sub_out_dir = fullfile(participant_dir, sub_id);
                if ~exist(sub_out_dir, 'dir'), mkdir(sub_out_dir); end

                % if exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_level1_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']),'file') && ...
                %         exist(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_level1_', hemisphere, '_', Methods{i_method}, '.mat']),'file')
                %     continue
                % end

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

                T = size(Ux, 3);

                % Calculate Divergence and Curl
                [DivMap, CurlMap] = velocity_derivatives(Ux, Uy, V_mask, 3, 1, 0);

                % --- Save voxel-wise Pattern density ---

                if cfg.do_pattern_voxel

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

                        SourceMap_Task = mean(SoM(:,:,Events(i_sub,1:T)==2),3,'omitnan');
                        SinkMap_Task = mean(SiM(:,:,Events(i_sub,1:T)==2),3,'omitnan');
                        CCVortexMap_Task = mean(CCM(:,:,Events(i_sub,1:T)==2),3,'omitnan');
                        CVortexMap_Task = mean(CM(:,:,Events(i_sub,1:T)==2),3,'omitnan');

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

                        Map = {SoM,SiM,CCM,CM};
                        for i_pattern = 1:numel(Map)
                            data = grid_to_vector(Map{i_pattern},V_mask)';
                            for i_compar = size(Comparison,1):-1:1
                                first_event = event_regressors.(Comparison{i_compar,1});
                                second_event = event_regressors.(Comparison{i_compar,2});
                                design_mat = [ones(T,1), first_event(i_sub,1:T)', second_event(i_sub,1:T)'];
                                cope(:,:,i_compar,i_pattern) = fit_gls_ar1(data, design_mat, c);
                            end
                        end
                        save(fullfile(sub_out_dir, [sub_id '_' Experiment '_numPattern_level1_' hemisphere '_' Methods{i_method} '.mat']),...
                            'cope');
                    end
                end

                % Calculate binarized Divergence and Curl
                [DivMap, CurlMap] = velocity_derivatives(Ux, Uy, V_mask, 3, 1, 1, 2);

                % Calculate Order Parameter
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
                if cfg.do_order_param_global
                    order_param_global = squeeze(mean(order_param, [1 2], 'omitnan'));
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_global', '-v7.3');
                end

                % Regional
                if cfg.do_order_param_region
                    order_param_region = zeros(num_regions, T);
                    for i_region = 1:num_regions
                        mask_rep = repmat(V_label == i_region, 1, 1, T);
                        tmp_param = order_param;
                        tmp_param(~mask_rep) = nan;
                        order_param_region(i_region, :) = squeeze(mean(tmp_param, [1 2], 'omitnan'));
                    end
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_region', '-v7.3');
                end

                % Voxel-wise
                if cfg.do_order_param_voxel
                    if contains(Experiment, 'rfMRI')
                        order_param_mean = mean(order_param, 3, 'omitnan');
                        save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_mean', '-v7.3');
                    else
                        if strcmp(Session{i_session},'MOTOR')
                            order_param_Fixation = mean(order_param(:, :, Events(i_sub, 1:T) == 0), 3, 'omitnan');
                            order_param_Task = mean(order_param(:, :, Events(i_sub, 1:T) == 1), 3, 'omitnan');
                            save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_Fixation', 'order_param_Task', '-v7.3');
                        else
                            order_param_Fixation = mean(order_param(:, :, Events(i_sub, 1:T) == 0), 3, 'omitnan');
                            order_param_Control = mean(order_param(:, :, Events(i_sub, 1:T) == 1), 3, 'omitnan');
                            order_param_Task = mean(order_param(:, :, Events(i_sub, 1:T) == 2), 3, 'omitnan');
                            save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']), 'order_param_Fixation', 'order_param_Control', 'order_param_Task', '-v7.3');
                        end

                        % effect size
                        data = grid_to_vector(order_param,V_mask)'; % T*V
                        cope = [];
                        for i_compar = size(Comparison,1):-1:1
                            first_event = event_regressors.(Comparison{i_compar,1});
                            second_event = event_regressors.(Comparison{i_compar,2});
                            design_mat = [ones(T,1), first_event(i_sub,1:T)', second_event(i_sub,1:T)'];
                            cope(:,:,i_compar) = fit_gls_ar1(data, design_mat, c);
                        end
                        save(fullfile(sub_out_dir, [sub_id '_' Experiment '_orderParam_level1_' hemisphere '_' Methods{i_method} '.mat']),...
                            'cope');
                    end
                end

                % --- 2 Save Pattern Sizes ---

                if cfg.do_pattern_size_global
                    SizeSource = zeros(1, T); SizeSink = zeros(1, T);
                    SizeCCVortex = zeros(1, T); SizeCVortex = zeros(1, T);
                    for t = 1:T
                        % Global size
                        [SizeSource(1,t), SizeSink(1,t)] = measure_motif_sizes(DivMap(:,:,t));
                        [SizeCCVortex(1,t), SizeCVortex(1,t)] = measure_motif_sizes(CurlMap(:,:,t));
                    end
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                        'SizeSource', 'SizeSink', 'SizeCCVortex', 'SizeCVortex', '-v7.3');
                end

                if cfg.do_pattern_size_region
                    SizeSource_region = zeros(num_regions, T); SizeSink_region = zeros(num_regions, T);
                    SizeCCVortex_region = zeros(num_regions, T); SizeCVortex_region = zeros(num_regions, T);
                    for t = 1:T
                        % Regional size
                        for i_region = 1:num_regions
                            region_mask = (V_label == i_region);

                            tmp_DivMap = DivMap(:,:,t); tmp_DivMap(~region_mask) = 0;
                            tmp_CurlMap = CurlMap(:,:,t); tmp_CurlMap(~region_mask) = 0;

                            [SizeSource_region(i_region, t), SizeSink_region(i_region, t)] = measure_motif_sizes(tmp_DivMap);
                            [SizeCCVortex_region(i_region, t), SizeCVortex_region(i_region, t)] = measure_motif_sizes(tmp_CurlMap);
                        end
                    end
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                        'SizeSource_region', 'SizeSink_region', 'SizeCCVortex_region', 'SizeCVortex_region', '-v7.3');
                end

                % --- 3 Save Pattern Quantity ---

                if cfg.do_pattern_quantity_global
                    NumSource = zeros(1, T); NumSink = zeros(1, T);
                    NumCCVortex = zeros(1, T); NumCVortex = zeros(1, T);
                    for t = 1:T
                        % Global size
                        [NumSource(1,t), NumSink(1,t)] = count_motifs(DivMap(:,:,t));
                        [NumCCVortex(1,t), NumCVortex(1,t)] = count_motifs(CurlMap(:,:,t));
                    end
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                        'NumSource', 'NumSink', 'NumCCVortex', 'NumCVortex', '-v7.3');
                end

                if cfg.do_pattern_quantity_region
                    NumSource_region = zeros(num_regions, T); NumSink_region = zeros(num_regions, T);
                    NumCCVortex_region = zeros(num_regions, T); NumCVortex_region = zeros(num_regions, T);
                    for t = 1:T
                        % Regional size
                        for i_region = 1:num_regions
                            region_mask = (V_label == i_region);

                            tmp_DivMap = DivMap(:,:,t); tmp_DivMap(~region_mask) = 0;
                            tmp_CurlMap = CurlMap(:,:,t); tmp_CurlMap(~region_mask) = 0;

                            [NumSource_region(i_region, t), NumSink_region(i_region, t)] = count_motifs(tmp_DivMap);
                            [NumCCVortex_region(i_region, t), NumCVortex_region(i_region, t)] = count_motifs(tmp_CurlMap);
                        end
                    end
                    save(fullfile(sub_out_dir, [sub_id, '_', Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']), ...
                        'NumSource_region', 'NumSink_region', 'NumCCVortex_region', 'NumCVortex_region', '-v7.3');
                end
            end
        end

    elseif strcmp(save_method, 'group')

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
                if cfg.do_order_param_voxel && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_order_param_voxel
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
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_level1_', hemisphere, '_', Methods{i_method}, '.mat']);
                if cfg.do_order_param_voxel && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_order_param_voxel
                    tmp = load(tmp_file);

                    if ~contains(Experiment,'rfMRI')
                        if i_sub==1
                            order_param_cope = zeros([num_subjects nnz(V_mask==1) size(Comparison,1)]);
                        end
                        order_param_cope(i_sub,:,:) = tmp.cope;
                    end
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

                tmp_file = fullfile(participant_dir, sub_id, [sub_id, '_', Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                if cfg.do_order_param_global && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_order_param_global
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
                if cfg.do_order_param_region && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_order_param_region
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
                if cfg.do_pattern_quantity_global && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_quantity_global
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
                if cfg.do_pattern_quantity_region && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_quantity_region
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
                if cfg.do_pattern_size_global && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_size_global
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
                if cfg.do_pattern_size_region && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_size_region
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
                if cfg.do_pattern_voxel && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_voxel
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
                if cfg.do_pattern_voxel && ~isfile(tmp_file)
                    error('fmri_task_metrics:MissingInput','Missing required participant file: %s',tmp_file);
                end
                if exist(tmp_file, 'file') && cfg.do_pattern_voxel
                    tmp = load(tmp_file);

                    if ~contains(Experiment,'rfMRI')
                        if i_sub==1
                            pattern_cope = zeros([num_subjects nnz(V_mask==1) size(Comparison,1) 4]);
                        end
                        pattern_cope(i_sub,:,:,:) = tmp.cope;
                    end
                else
                    warning(['File missing for Subject ', sub_id, ': ', tmp_file]);
                end

            end

            % 保存该 Method 和 Task 下的组级矩阵
            if cfg.do_order_param_global
                group_save_name = fullfile(group_dir, [Experiment, '_orderParam_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'order_param_global', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_order_param_region
                group_save_name = fullfile(group_dir, [Experiment, '_orderParam_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'order_param_region', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_order_param_voxel
                if contains(Experiment,'rfMRI')
                    group_save_name = fullfile(group_dir, [Experiment, '_orderParam_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
                    save(group_save_name, 'order_param', '-v7.3');
                else
                    group_save_name = fullfile(group_dir, [Experiment, '_orderParam_voxel_ttest_', hemisphere, '_', Methods{i_method}, '.mat']);
                    % voxel-wise paired t-test
                    t_values = [];
                    p_values = [];
                    for i_compar = 1:size(Comparison,1)
                        Map_1 = select_condition_map(Comparison{i_compar,1}, order_param_Fixation, order_param_Control, order_param_Task);
                        Map_2 = select_condition_map(Comparison{i_compar,2}, order_param_Fixation, order_param_Control, order_param_Task);
                        Map_1 = grid_to_vector(Map_1, V_mask);
                        Map_2 = grid_to_vector(Map_2, V_mask);
                        for i_vertex = 1:size(Map_1,1)
                            [~, p_values(i_vertex, i_compar), ~, stats] = ttest(Map_1(i_vertex, :)', Map_2(i_vertex, :)');
                            t_values(i_vertex, i_compar) = stats.tstat;
                        end
                    end
                    save(group_save_name, 't_values', 'p_values', '-v7.3');
                end
                disp(['Saved: ', group_save_name]);

                group_save_name = fullfile(group_dir, [Experiment, '_orderParam_cohensd_', hemisphere, '_', Methods{i_method}, '.mat']);
                if ~contains(Experiment,'rfMRI')
                    % cope: nSub × V × nCompar
                    group_mean = squeeze(mean(order_param_cope, 1, 'omitnan'));          % V × nCompar
                    group_sd   = squeeze(std(order_param_cope, 0, 1, 'omitnan'));        % V × nCompar, sample SD
                    cohens_d   = group_mean ./ (group_sd + eps);  % V × nCompar x 4
                    save(group_save_name, 'cohens_d', 'Comparison', '-v7.3');
                end
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_pattern_quantity_global
                group_save_name = fullfile(group_dir, [Experiment, '_numPattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'NumSource', 'NumSink', 'NumCCVortex', 'NumCVortex', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_pattern_quantity_region
                group_save_name = fullfile(group_dir, [Experiment, '_numPattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'NumSource_region', 'NumSink_region', 'NumCCVortex_region', 'NumCVortex_region', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_pattern_size_global
                group_save_name = fullfile(group_dir, [Experiment, '_sizePattern_global_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'SizeSource', 'SizeSink', 'SizeCCVortex', 'SizeCVortex', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_pattern_size_region
                group_save_name = fullfile(group_dir, [Experiment, '_sizePattern_region-', parcellation, '_', hemisphere, '_', Methods{i_method}, '.mat']);
                save(group_save_name, 'SizeSource_region', 'SizeSink_region', 'SizeCCVortex_region', 'SizeCVortex_region', '-v7.3');
                disp(['Saved: ', group_save_name]);
            end

            if cfg.do_pattern_voxel
                if contains(Experiment,'rfMRI')
                    group_save_name = fullfile(group_dir, [Experiment, '_numPattern_voxel_', hemisphere, '_', Methods{i_method}, '.mat']);
                    save(group_save_name, 'SourceMap', 'SinkMap', 'CCVortexMap', 'CVortexMap', '-v7.3');
                else
                    group_save_name = fullfile(group_dir, [Experiment, '_numPattern_voxel_ttest_', hemisphere, '_', Methods{i_method}, '.mat']);
                    % voxel-wise paired t-test
                    t_values = [];
                    p_values = [];
                    for i_compar = 1:size(Comparison,1)
                        for i_pattern = 4:-1:1
                            pattern_maps = {SourceMap_Fixation, SourceMap_Control, SourceMap_Task;
                                SinkMap_Fixation, SinkMap_Control, SinkMap_Task;
                                CCVortexMap_Fixation, CCVortexMap_Control, CCVortexMap_Task;
                                CVortexMap_Fixation, CVortexMap_Control, CVortexMap_Task};
                            names = {'Fixation','Control','Task'};
                            Map_1 = pattern_maps{i_pattern,strcmp(names,Comparison{i_compar,1})};
                            Map_2 = pattern_maps{i_pattern,strcmp(names,Comparison{i_compar,2})};
                            Map_1 = grid_to_vector(Map_1, V_mask);
                            Map_2 = grid_to_vector(Map_2, V_mask);
                            for i_vertex = 1:size(Map_1,1)
                                [~, p_values(i_vertex, i_pattern, i_compar), ~, stats] = ttest(Map_1(i_vertex, :)', Map_2(i_vertex, :)');
                                t_values(i_vertex, i_pattern, i_compar) = stats.tstat;
                            end
                        end
                    end
                    save(group_save_name, 't_values', 'p_values', '-v7.3');
                end
                disp(['Saved: ', group_save_name]);

                group_save_name = fullfile(group_dir, [Experiment, '_numPattern_cohensd_', hemisphere, '_', Methods{i_method}, '.mat']);
                if ~contains(Experiment,'rfMRI')
                    % cope: nSub × V × nCompar x 4
                    group_mean = squeeze(mean(pattern_cope, 1, 'omitnan'));          % V × nCompar x 4
                    group_sd   = squeeze(std(pattern_cope, 0, 1, 'omitnan'));        % V × nCompar x 4, sample SD
                    cohens_d   = group_mean ./ (group_sd + eps);  % V × nCompar x 4
                    save(group_save_name, 'cohens_d', 'Comparison', '-v7.3');
                end
                disp(['Saved: ', group_save_name])
            end

        end
        disp('Group level aggregation completed!');
    end

end

end
function value = select_condition_map(name,fixation,control,task)
switch name
    case 'Fixation', value=fixation;
    case 'Control', value=control;
    case 'Task', value=task;
    otherwise, error('FMRI:Condition','Unknown condition: %s.',name);
end
end
