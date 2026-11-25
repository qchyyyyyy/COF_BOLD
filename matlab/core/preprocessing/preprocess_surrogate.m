function is_valid = preprocess_surrogate(nii_file, proc_dir,surrogateMethod)
% preprocess_rest: map fmri-time-series in cifti format into meshgrid with cord_file
% and generate preprocessed-mesh-grid data in proc_dir

% Input:
% nii_file: single person fmri-time-series-data; for each time point, two-dimensional cortex surface
% proc_dir: directory where generated mesh-grid data stored in

% Output:
% data.BrainImg_lh: preprocessed mesh-grid left fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))
% data.BrainImg_rh: preprocessed mesh-grid right fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))

spmRoot=getenv('SPM_ROOT');
assert(~isempty(spmRoot)&&exist(spmRoot,'dir')==7,'Set SPM_ROOT to the SPM installation.');
addpath(spmRoot);

load(template_path('coordinate_cii_v2.mat'), 'x_left', 'y_left', 'x_right', 'y_right', 'cii_mask_left', 'cii_mask_right', 'i_invalid_left', 'i_invalid_right');

% if exist(nii_file,'file')==0 || exist(fullfile(proc_dir,'TC_mesh_r_f_zscore.mat'),'file')~=0
%     return
% end

is_valid = true;
if ~exist(nii_file,'file')
    disp([nii_file ' not exist']);
    is_valid = false;
    return
end

if exist([proc_dir,'TC_mesh_r_f_zscore.mat'],'file')
    disp([proc_dir,'TC_mesh_r_f_zscore.mat already exist']);
    return
end

try
    fMRI2d = ft_read_cifti(nii_file);
catch ME
    disp([nii_file ' ft_read_cifti error:' ME.message]);
    is_valid = false;
    return;
end


%% preparations
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
downSRate = 2 ;
x_target_cood = -250:downSRate:250;
y_target_cood = -200:downSRate:200;
[xi,yi] = meshgrid(x_target_cood,y_target_cood);

%% joint processing of bandpass filter and global signal regression
fsTem = 1/TR;
fLow = 0.01 ;
fHigh = 0.1 ;
idx_left_end = size(TC_lh,1);

Y = [TC_lh' TC_rh'];clearvars TC_lh TC_rh

Y = Y - mean(Y);

GS = mean(Y,2,'omitnan');

M = [Y , GS];
M_bp = bandpass(M,[fLow fHigh],fsTem);clearvars M
Y_bp = M_bp(:,1:size(Y,2));
X_bp = M_bp(:,size(Y,2)+1:end);clearvars Y M_bp GS

beta = X_bp \ Y_bp;
Y = Y_bp - X_bp*beta; clearvars Y_bp X_bp beta

T = size(Y,1);

if TR>0.8
    Y = interpft(Y,round(TR/0.72*T),1);
end

T = size(Y,1);
%% interpolate data onto 2-D grid

TC_lh = Y(:,1:idx_left_end)';
TC_rh = Y(:,idx_left_end+1:end)';clearvars Y
TC_lh_mesh = zeros(length(y_target_cood), length(x_target_cood), T);
TC_rh_mesh = TC_lh_mesh;
% Restore warning state when this function returns, including on failure.
warningState=warning; warningCleanup=onCleanup(@() warning(warningState)); %#ok<NASGU>
warning('off');

for i_timepoint = 1:T
    % left brain
    F_TC_lh = scatteredInterpolant(x_left(~isnan(TC_lh(:,1))),y_left(~isnan(TC_lh(:,1))),TC_lh(~isnan(TC_lh(:,1)),i_timepoint));
    TC_lh_mesh(:,:,i_timepoint) = F_TC_lh(xi, yi).*cii_mask_left;
    % right brain
    F_TC_rh = scatteredInterpolant(x_right(~isnan(TC_rh(:,1))),y_right(~isnan(TC_rh(:,1))),TC_rh(~isnan(TC_rh(:,1)),i_timepoint));
    TC_rh_mesh(:,:,i_timepoint) = F_TC_rh(xi, yi).*cii_mask_right;
end
clearvars TC_lh TC_rh

TC_lh_mesh = zscore(TC_lh_mesh,[],3);
TC_rh_mesh = zscore(TC_rh_mesh,[],3);

%% generate surrogate data
TC = cat(2,TC_lh_mesh,TC_rh_mesh);
mask = [cii_mask_left cii_mask_right];
clearvars TC_lh TC_rh TC_lh_mesh TC_rh_mesh F_TC_lh F_TC_rh
TC(:,251,:)=[];
mask(:,251)=[];

switch surrogateMethod
    case '3D'
        TC = fft_3d(TC);
    case '2D'
        TC = fft_2d(grid_to_vector(TC,mask));
        TC = vector_to_grid(TC,mask);
end

TC_lh_mesh = TC(:,1:251,:);
TC_rh_mesh = TC(:,251:end,:);


%% save data
data.BrainImg_lh = TC_lh_mesh;
data.BrainImg_rh = TC_rh_mesh;
data.F_TR = TR;
data.F_band = [fLow fHigh];
data.Normalization = 'zscore';%'none' 'demean'
data.R_globalSig = 'yes';
data.R_deconv = 'no';

if ~exist(proc_dir,'dir')
    mkdir(proc_dir);
end

save([proc_dir,'TC_surrogate_' surrogateMethod '.mat'], 'data');
disp([proc_dir,' BrainImg is sucessfully preprocessed and saved']);

end

function dataOut = fft_3d(dataIn)

midPoint = floor(size(dataIn)/2)+1;
indMid = sub2ind(size(dataIn),midPoint(1),midPoint(2),midPoint(3));
freqData = fftshift(fftn(dataIn)); % 3D fourier transform
phaseData = angle(freqData);
absData = abs(freqData);
phaseDataRand = phaseData;
phaseDataRand(1:indMid-1) = rand(indMid-1,1).*2*pi-pi; % random phase values -pi~pi
for iInd = 1:indMid-2
    phaseDataRand(indMid+iInd) = -phaseDataRand(indMid-iInd) ;
end
freqDataNew = ifftshift(absData.*exp(1i.*phaseData + 1i*phaseDataRand)); % random phase shuflling
dataOut = real(ifftn(freqDataNew)) ; % 3D inverse fourier tranform

end

function dataOut = fft_2d(dataIn)

midPoint = floor(size(dataIn,2)/2)+1;

freqData = fftshift(fft(dataIn,[],2)); % 2D fourier transform
phaseData = angle(freqData);
absData = abs(freqData);
phaseDataRand = phaseData;
phaseDataRand(:,1:midPoint-1) = rand(size(phaseData,1),midPoint-1).*2*pi-pi; % random phase values -pi~pi
for iInd = 1:midPoint-2
    phaseDataRand(:,midPoint+iInd) = -phaseDataRand(:,midPoint-iInd) ;
end
freqDataNew = ifftshift(absData.*exp(1i.*phaseData + 1i*phaseDataRand)); % random phase shuflling
dataOut = real(ifft(freqDataNew,[],2)) ; % 2D inverse fourier tranform

end
