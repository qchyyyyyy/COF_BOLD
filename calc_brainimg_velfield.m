clear;clc
addpath(genpath('functions'))

dataPath = 'data/fMRI/rfMRI_REST1_LR/';
% dataPath = 'data/fMRI/tfMRI_WM_LR/';
% dataPath = 'data/fMRI/tfMRI_WM_RL/';
% dataPath = 'data/fMRI/tfMRI_RELATIONAL_LR/';
% dataPath = 'data/fMRI/tfMRI_RELATIONAL_RL/';
% dataPath = 'data/fMRI/tfMRI_LANGUAGE_LR/';
% dataPath = 'data/fMRI/tfMRI_MOTOR_LR/';

File = dir(fullfile(dataPath));
SubjectList = cell(size(File));
for m = 1:length(File)
   SubjectList{m}=File(m).name;
end
SubjectList(1:3)=[];
%% data preprocessing
for sub = 1:length(SubjectList)
    disp(['sub = ' num2str(sub)]);
    preproc_2D([dataPath,SubjectList{sub},'/rfMRI_REST1_LR_Atlas_hp2000_clean.dtseries.nii'],...
        [dataPath,SubjectList{sub},'/']);
end

%% VVFs estimation via COF
load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_mask_right');
tau = 5e-4;
rho = 1e-8;
sigma = 8;% 16mm

for sub = length(SubjectList):-1:1
    fMRI{sub}.path = [dataPath,SubjectList{sub}, '/'];
    fMRI{sub}.fileName = 'TC_mesh_r_f_zscore.mat';
    fMRI{sub}.hemisphere = 'lh';
end

for sub = 1:length(SubjectList)
    disp(['sub = ' num2str(sub)]);
    extract_velocity_field_hilbert(fMRI{sub},cii_mask_left,tau,rho,sigma);
end
%%
for sub = length(SubjectList):-1:1
    fMRI{sub}.hemisphere = 'rh';
end

for sub = 1:length(SubjectList)
    disp(['sub = ' num2str(sub)]);
    extract_velocity_field_hilbert(fMRI{sub},cii_mask_right,tau,rho,sigma);
end