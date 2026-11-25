function is_valid = preprocess_task(nii_file, proc_dir,TR,cii_mask_left,cii_mask_right,x_target_cood,y_target_cood,x_left,y_left,x_right,y_right)
% PREPROC_2D_TFMRI: map fmri-time-series in cifti format into meshgrid with cord_file
% and generate preprocessed-mesh-grid data in proc_dir

% Input:
% nii_file: single person fmri-time-series-data; for each time point, two-dimensional cortex surface
% proc_dir: directory where generated mesh-grid data stored in

% Output:
% data.BrainImg_lh: preprocessed mesh-grid left fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))
% data.BrainImg_rh: preprocessed mesh-grid right fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))

if ~exist('TR','var')
    TR = 0.72;
end

if ~exist('x_target_cood','var')
    downSRate = 2 ;
    x_target_cood = -250:downSRate:250;
    y_target_cood = -200:downSRate:200;
end

if ~exist('cii_mask_left','var')
    load(template_path('coordinate_cii_v2.mat'), 'cii_mask_left', 'cii_mask_right');
end

if ~exist('x_left','var')
    surf_left = gifti(template_path('HCP_S1200_Atlas_Z4_pkXDZ/S1200.L.flat.32k_fs_LR.surf.gii'));
    surf_right = gifti(template_path('HCP_S1200_Atlas_Z4_pkXDZ/S1200.R.flat.32k_fs_LR.surf.gii'));
    x_left = double(surf_left.vertices(:,1));
    y_left = double(surf_left.vertices(:,2));
    x_right = double(surf_right.vertices(:,1));
    y_right = double(surf_right.vertices(:,2));
end

is_valid = true;
if ~exist(nii_file,'file')
    disp([nii_file ' not exist']);
    is_valid = false;
    return
end

preprocessed_file = fullfile(proc_dir, 'TC_mesh_r_f_zscore.mat');
if exist(preprocessed_file,'file')
    disp([preprocessed_file ' already exists']);
    return
end

try
    fMRI2d = ft_read_cifti(nii_file);
catch ME
    disp([nii_file ' ft_read_cifti error:' ME.message]);
    is_valid = false;
    return;
end

% preparations
TC_2D = fMRI2d.dtseries;
TC_2D(sum(TC_2D,2)==0,:)=nan;
TC_lh = TC_2D(fMRI2d.brainstructure==1,:);
TC_rh = TC_2D(fMRI2d.brainstructure==2,:);
i_invalid_left = isnan(TC_lh(:,1));
i_invalid_right = isnan(TC_rh(:,1));
TC_lh(i_invalid_left,:) = [];
TC_rh(i_invalid_right,:) = [];
x_left(i_invalid_left,:) = [];
y_left(i_invalid_left,:) = [];
x_right(i_invalid_right,:) = [];
y_right(i_invalid_right,:) = [];
clear TC_2D
[xi,yi] = meshgrid(x_target_cood,y_target_cood);

% high-pass filter and global signal regression
idx_left_end = size(TC_lh,1);

Y = [TC_lh' TC_rh']; clearvars TC_lh TC_rh

Y = Y - mean(Y);

GS = mean(Y,2,'omitnan');

T = size(Y,1);

K0 = struct('RT',TR,'HParam',128,'row',1:T);
K = spm_filter(K0);
Y_hp = spm_filter(K, Y);
X_hp = spm_filter(K, GS);

beta = (X_hp' * X_hp) \ (X_hp' * Y_hp);
Y  = Y_hp - X_hp * beta;

clearvars GS beta Y_hp X_hp

h   = spm_hrf(TR);
H   = fft(h, size(Y,1));
lam = 0.2;
W   = conj(H)./(abs(H).^2+lam);
Y = real(ifft(fft(Y).*W));

fHigh = min(0.25,1./TR./2);
Y = lowpass(Y, fHigh, 1/TR);

if TR>0.8
    Y = interpft(Y,round(TR/0.72*T),1);
end

T = size(Y,1);

TC_lh = Y(:,1:idx_left_end)';
TC_rh = Y(:,idx_left_end+1:end)'; clearvars Y


% interpolate data onto 2-D grid
TC_lh_mesh = zeros(length(y_target_cood), length(x_target_cood), T);
TC_rh_mesh = TC_lh_mesh;
T = size(TC_lh,2);
% Restore warning state when this function returns, including on failure.
warningState=warning; warningCleanup=onCleanup(@() warning(warningState)); %#ok<NASGU>
warning('off');

for i_timepoint = 1:T
    % left brain
    F_TC_lh = scatteredInterpolant(x_left,y_left,TC_lh(:,i_timepoint));
    TC_lh_mesh(:,:,i_timepoint) = F_TC_lh(xi, yi).*cii_mask_left;
    % right brain
    F_TC_rh = scatteredInterpolant(x_right,y_right,TC_rh(:,i_timepoint));
    TC_rh_mesh(:,:,i_timepoint) = F_TC_rh(xi, yi).*cii_mask_right;
end
clearvars TC_lh TC_rh

TC_lh_mesh = zscore(TC_lh_mesh,[],3);
TC_rh_mesh = zscore(TC_rh_mesh,[],3);

% save data
data.BrainImg_lh = TC_lh_mesh;
data.BrainImg_rh = TC_rh_mesh;
data.F_TR = TR;
data.F_band = [1/128 fHigh];
data.Normalization = 'zscore';%'none' 'demean'
data.R_globalSig = 'yes';

if ~exist(proc_dir,'dir')
    mkdir(proc_dir);
end
save(preprocessed_file, 'data');
disp([proc_dir,' BrainImg is sucessfully preprocessed and saved']);

end
