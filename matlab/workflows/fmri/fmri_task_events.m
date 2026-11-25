function fmri_task_events(taskIndices)
% Task indices: 1 emotion, 2 gambling, 3 language, 5 relational, 6 social, 7 WM.
if nargin<1, taskIndices=[1 2 3 5 6 7]; end
assert(~isempty(getenv('COF_DATA_ROOT')),'Set COF_DATA_ROOT to external derived data.');

% Run setup('fmri') before this event exporter.

Site = 'HCP-YA';
Type = 'tfMRI';
PE = 'LR';

raw_data_dir = getenv('HCP_TASK_RAW_ROOT');
assert(~isempty(raw_data_dir)&&exist(raw_data_dir,'dir')==7,'Set HCP_TASK_RAW_ROOT.');

Session = {'EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM'};
Session_name = {'Emotion', 'Gambling', 'Language', 'Motor', 'Relational', 'Social', 'Working_Memory'};

for i_session = taskIndices
    
    raw_data_dir_task = fullfile(raw_data_dir, Session_name{i_session}, 'preproc');
    Experiment = [Type '_' Session{i_session} '_' PE];
    data_dir = fullfile(getenv('COF_DATA_ROOT'), Site, Experiment);
    save_dir = fullfile(getenv('COF_DATA_ROOT'), 'results', Site, Experiment, 'group'); [~, ~] = mkdir(save_dir);

    subjectList = dir(fullfile(raw_data_dir_task));
    subjectList = {subjectList.name};
    subjectList(1:2) = [];

    for i_sub = numel(subjectList):-1:1
        img_file = fullfile(raw_data_dir_task, subjectList{i_sub}, 'MNINonLinear', 'Results', Experiment, [Experiment '_Atlas_MSMAll.dtseries.nii']);
        if ~exist(img_file,'file')
            subjectList(i_sub) = [];
        end
    end
    save(fullfile(save_dir, [Experiment '_SubList.mat']),'subjectList');
    num_subjects = numel(subjectList);

    load(fullfile(data_dir,subjectList{1},'VelField_lh_COF.mat'),'VelField');
    T = size(VelField.Ux,3); clearvars VelField
    TR = 0.72;

    EVs_file = cell(num_subjects, 1);
    Events = zeros(num_subjects, T);
    Stim_Type_subs = [];

    for i_sub = num_subjects:-1:1
        disp(['sub ' num2str(i_sub)]);

        sub_id = subjectList{i_sub};
        EVs_dir = fullfile(raw_data_dir_task, sub_id, 'MNINonLinear', 'Results', Experiment, 'EVs');

        if ~exist(EVs_dir,'dir')
            continue
        end
        EVs_file_sub = table('Size', [0 4], ...
            'VariableTypes', {'string', 'double', 'double', 'double'}, ...
            'VariableNames', {'Conditions', 'Onset', 'End', 'Amplitude'});

        Session = {'EMOTION', 'GAMBLING', 'LANGUAGE', 'MOTOR', 'RELATIONAL', 'SOCIAL', 'WM'};


        switch Session{i_session}
            case 'EMOTION'
                condition_names = {'fear', 'neut'};
            case 'GAMBLING'
                condition_names = {'win', 'loss'};
            case 'LANGUAGE'
                condition_names = {'story', 'math'};
            case 'MOTOR'
                condition_names = {'lf', 'lh', 'rf', 'rh', 't'};
            case 'RELATIONAL'
                condition_names = {'relation', 'match'};
            case 'SOCIAL'
                condition_names = {'mental', 'rnd'};
            case 'WM'
                condition_names = {'0bk_body', '0bk_faces', '0bk_places', '0bk_tools', '2bk_body', '2bk_faces', '2bk_places', '2bk_tools'};
        end

        for i_condition = 1:numel(condition_names)
            try
                tmp_EVs_file = readtable(fullfile(EVs_dir, [condition_names{i_condition} '.txt']), 'ReadVariableNames', false, 'Delimiter', '\t');
            catch
                continue
            end
            tmp_EVs_file
            tmp_EVs_file.Conditions = repmat(condition_names{i_condition}, size(tmp_EVs_file, 1), 1);
            tmp_EVs_file.Properties.VariableNames = {'Onset', 'End', 'Amplitude', 'Conditions'};
            tmp_EVs_file.End = tmp_EVs_file.End + tmp_EVs_file.Onset;

            EVs_file_sub = [EVs_file_sub; tmp_EVs_file(:, [4 1 2 3])];
        end

        [~, idx] = sort(EVs_file_sub.Onset);
        EVs_file_sub = EVs_file_sub(idx, :);

        EVs_file_sub.Onset = EVs_file_sub.Onset./TR;
        EVs_file_sub.End = EVs_file_sub.End./TR;

        EVs_file{i_sub} = EVs_file_sub;


        n_blocks = length(EVs_file_sub.Onset);
        for i_block = 1:n_blocks
            startTime = max(1,round(EVs_file_sub.Onset(i_block)));
            endTime = min(T,round(EVs_file_sub.End(i_block)));
            if strcmp(Session{i_session},'MOTOR')
                label_condition = 1;
            elseif strcmp(Session{i_session},'WM')
                label_condition = contains(EVs_file_sub.Conditions(idx(i_block)), '2bk') + 1;
            else
                label_condition = contains(EVs_file_sub.Conditions(idx(i_block)), condition_names{1}) + 1;
            end
            Stim_Type_subs(i_block,:,i_sub) = [startTime endTime label_condition];
            Events(i_sub, startTime:endTime) = label_condition;

        end
    end

    save(fullfile(save_dir, [Experiment '_EVs.mat']), 'Events', 'EVs_file','Stim_Type_subs','condition_names');

end

end
