function fmri_gambling_events()
% function fmri_gambling_events(Site,Experiment,data_dir,subjectList,TR)
clear; clc;
assert(~isempty(getenv('COF_DATA_ROOT')),'Set COF_DATA_ROOT to external derived data.');

Site = 'HCP_3T';
Experiment = 'tfMRI_GAMBLING_LR';
raw_data_dir = getenv('HCP_GAMBLING_PREPROC_ROOT');
assert(~isempty(raw_data_dir)&&exist(raw_data_dir,'dir')==7,'Set HCP_GAMBLING_PREPROC_ROOT.');
data_dir = fullfile(getenv('COF_DATA_ROOT'), Site, Experiment);
subjectList = dir(fullfile(data_dir));
subjectList = {subjectList.name};
subjectList(1:2) = [];

TR = 0.72;
T = 252;
results_dir = fullfile(getenv('COF_DATA_ROOT'), 'results', Site, Experiment, 'group');
[~, ~] = mkdir(results_dir);
% Run setup('fmri') before this event exporter.
num_subjects = numel(subjectList);

EVs_file = cell(num_subjects, 1);
Events = zeros(num_subjects, T);

for i_sub = num_subjects:-1:1
    disp(['sub ' num2str(i_sub)]);

    sub_id = subjectList{i_sub};
    EVs_dir = fullfile(raw_data_dir, sub_id, 'MNINonLinear', 'Results', Experiment, 'EVs');
    
    condition_names = {'win.txt','loss.txt'};
    
    EVs_file_sub = table('Size', [0 4], ...
              'VariableTypes', {'string', 'double', 'double', 'double'}, ...
              'VariableNames', {'Conditions', 'Onset', 'End', 'Amplitude'});
    for i_condition = 1:numel(condition_names)

        tmp_EVs_file = readtable(fullfile(EVs_dir, condition_names{i_condition}));
        tmp_name = split(condition_names{i_condition}, '.');
        tmp_EVs_file.Conditions = repmat(tmp_name{1}, size(tmp_EVs_file, 1), 1);
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
        startTime = EVs_file_sub.Onset(i_block);
        endTime = EVs_file_sub.End(i_block);

        label_condition = contains(EVs_file_sub.Conditions(idx(i_block)), 'win') + 1;
        Events(i_sub, max(1, round(startTime)):round(endTime)) = label_condition;
        Stim_Type_subs(i_block,:,i_sub) = [startTime endTime label_condition];
    end
end

save(fullfile(results_dir, [Experiment '_EVs.mat']), 'Events', 'EVs_file','Stim_Type_subs');

end
