function fmri_bold_copes(save_method, i_batch, batch_size)

assert(~isempty(getenv('COF_DATA_ROOT')) && ~isempty(getenv('COF_RESULTS_ROOT')) && ~isempty(getenv('COF_TASK_OUTPUT_ROOT')),'Set external data, metadata and fresh output roots.');
assert(~strcmp(getenv('COF_RESULTS_ROOT'),getenv('COF_TASK_OUTPUT_ROOT')),'Input and output roots must differ.');

parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';

[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = load_cortical_atlas(parcellation,'L');
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = load_cortical_atlas(parcellation,'R');

mask = [mask_L(:,1:end-1) mask_R];
label = [label_L(:,1:end-1) label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
V_label_L = max(label_L(1:end-1,1:end-1),max(label_L(2:end,1:end-1),max(label_L(1:end-1,2:end),label_L(2:end,2:end))));
V_label_L(V_mask_L==0)=0;
V_label_R = max(label_R(1:end-1,1:end-1),max(label_R(2:end,1:end-1),max(label_R(1:end-1,2:end),label_R(2:end,2:end))));
V_label_R(V_mask_R==0)=0;
V_label = [V_label_L V_label_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R V_label_L V_label_R


Site = 'HCP-YA';
Type = 'tfMRI';
PE = 'LR';

raw_data_dir = getenv('HCP_TASK_RAW_ROOT'); % legacy unused motion-file construction
Session = {'EMOTION', 'GAMBLING', 'LANGUAGE', 'RELATIONAL', 'SOCIAL', 'WM'}; % , 'MOTOR'
Session_name = {'Emotion', 'Gambling', 'Language', 'Motor', 'Relational', 'Social', 'Working_Memory'};

for i_session = 1:numel(Session)

    
    Experiment = [Type '_' Session{i_session} '_' PE];
    
    disp(Experiment)
    
    raw_data_dir_task = fullfile(raw_data_dir, Session_name{i_session}, 'preproc');
    data_dir = fullfile(getenv('COF_DATA_ROOT'), Site, Experiment);
    participant_dir = fullfile(getenv('COF_TASK_OUTPUT_ROOT'),Site,Experiment,'participant');
    group_dir = fullfile(getenv('COF_TASK_OUTPUT_ROOT'),Site,Experiment,'group');
    metadata_dir=fullfile(getenv('COF_RESULTS_ROOT'),Site,Experiment,'group');
    if ~isfolder(group_dir),mkdir(group_dir);end

    load(fullfile(metadata_dir, [Experiment '_EVs.mat']),'Events');

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

    % Access the list of subject files
    load(fullfile(metadata_dir, [Experiment '_SubList.mat']), 'subjectList')
    num_subjects = length(subjectList);

    if strcmp(save_method, 'participant')

        if nargin >= 3 && ~isempty(i_batch) && ~isempty(batch_size)
            idx_sub = i_batch:batch_size:num_subjects;
        else
            idx_sub = 1:num_subjects;
        end

        for i_sub = idx_sub
            sub_id = subjectList{i_sub};
            disp(['sub ' num2str(i_sub)]);

            load(fullfile(data_dir, sub_id, 'TC_mesh_r_f_zscore.mat'),'data');
            data = cat(2,data.BrainImg_lh(:,1:end-1,:),data.BrainImg_rh);
            data = convn(data,ones(2,2,2)/8,'valid');
            data = grid_to_vector(data,V_mask)'; % T*V

            T = size(data,1);

            motion_file = fullfile(raw_data_dir_task, sub_id, 'MNINonLinear', 'Results', Experiment, 'EVs');

            c = [0; 1; -1]; % 3*1

            cope = []; varcope = []; tstat = []; beta_gls = [];
            for i_compar = size(Comparison,1):-1:1
                first_event = event_regressors.(Comparison{i_compar,1});
                second_event = event_regressors.(Comparison{i_compar,2});
                design_mat = [ones(T,1), first_event(i_sub,1:T)', second_event(i_sub,1:T)'];
                [cope(:,:,i_compar),varcope(:,:,i_compar),tstat(:,:,i_compar),beta_gls{i_compar}] = fit_gls_ar1(data, design_mat, c);
            end
            if ~exist(fullfile(participant_dir, sub_id),'dir'); mkdir(fullfile(participant_dir,sub_id)); end
            save(fullfile(participant_dir, sub_id, [sub_id '_' Experiment '_brainImg_level1_' hemisphere '.mat']),...
                'cope','varcope','tstat','beta_gls');

        end

    else
        cope = []; 
        % varcope = []; tstat = []; beta_gls = [];
        for i_sub = num_subjects:-1:1

            sub_id = subjectList{i_sub};
            
            disp(['sub ' num2str(i_sub)]);
            
            tmp = load(fullfile(participant_dir, sub_id, [sub_id '_' Experiment '_brainImg_level1_' hemisphere '.mat']),...
                'cope','varcope','tstat','beta_gls');
            cope(i_sub,:,:) = tmp.cope;
            % varcope(i_sub,:,:) = tmp.varcope;
            % tstat(i_sub,:,:) = tmp.tstat;
            % beta_gls(i_sub,:) = tmp.beta_gls;
        end

        % cope: nSub × V × nCompar
        group_mean = squeeze(mean(cope, 1, 'omitnan'));          % V × nCompar
        group_sd   = squeeze(std(cope, 0, 1, 'omitnan'));        % V × nCompar, sample SD
        cohens_d   = group_mean ./ (group_sd + eps);  % V × nCompar

        save(fullfile(group_dir, [Experiment '_level2_cohensd_' hemisphere '.mat']), ...
            'cohens_d', ...
            'Comparison', '-v7.3');
    end

end

end
