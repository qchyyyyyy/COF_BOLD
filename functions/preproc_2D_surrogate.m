function preproc_2D_surrogate(nii_file, proc_dir,surrogateMethod)
% preproc_2D: map fmri-time-series in cifti format into meshgrid with cord_file
% and generate preprocessed-mesh-grid data in proc_dir

% Input:
% nii_file: single person fmri-time-series-data; for each time point, two-dimensional cortex surface
% proc_dir: directory where generated mesh-grid data stored in

% Output:
% data.BrainImg_lh: preprocessed mesh-grid left fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))
% data.BrainImg_rh: preprocessed mesh-grid right fmri time series data (shape (length(y_target_cood), length(x_target_cood), timepoint))

load('data/templates/coordinate_cii_v2.mat', 'x_left', 'y_left', 'x_right', 'y_right', 'cii_mask_left', 'cii_mask_right', 'i_invalid_left', 'i_invalid_right');

% if exist(nii_file,'file')==0 || exist(fullfile(proc_dir,'TC_mesh_r_f_zscore.mat'),'file')~=0
%     return
% end

fMRI2d = ft_read_cifti(nii_file);

%% preparations
TC_2D = fMRI2d.dtseries;
TC_lh = TC_2D(fMRI2d.brainstructure==1,:);
TC_rh = TC_2D(fMRI2d.brainstructure==2,:);
TC_lh(i_invalid_left,:) = [];
TC_rh(i_invalid_right,:) = [];
clear TC_2D
downSRate = 2 ;
x_target_cood = -250:downSRate:250;
y_target_cood = -200:downSRate:200;
[xi,yi] = meshgrid(x_target_cood,y_target_cood);

%% interpolate data onto 2-D grid
T = size(TC_lh,2);
TC_lh_mesh = zeros(length(y_target_cood), length(x_target_cood), T-1);
TC_rh_mesh = TC_lh_mesh;
warning('off');

for i_timepoint = 1:T-1
    % left brain
    F_TC_lh = scatteredInterpolant(x_left(~isnan(TC_lh(:,1))),y_left(~isnan(TC_lh(:,1))),TC_lh(~isnan(TC_lh(:,1)),i_timepoint));
    TC_lh_mesh(:,:,i_timepoint) = F_TC_lh(xi, yi).*cii_mask_left;
    % right brain
    F_TC_rh = scatteredInterpolant(x_right(~isnan(TC_rh(:,1))),y_right(~isnan(TC_rh(:,1))),TC_rh(~isnan(TC_rh(:,1)),i_timepoint));
    TC_rh_mesh(:,:,i_timepoint) = F_TC_rh(xi, yi).*cii_mask_right;
end
clearvars TC_lh TC_rh
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
        TC = fft_2d(map2vec(TC,mask));
        TC = vec2map(TC,mask);
end

TC_lh_mesh = TC(:,1:251,:);
TC_rh_mesh = TC(:,251:end,:);
%% global signal regression
TC_lh_mesh = map2vec(TC_lh_mesh,cii_mask_left);
TC_rh_mesh = map2vec(TC_rh_mesh,cii_mask_right);

GS = mean(TC_lh_mesh,1,'omitnan');
beta = (TC_lh_mesh*GS')/(GS*GS');
TC_lh_mesh = TC_lh_mesh - beta*GS;

GS = mean(TC_rh_mesh,1,'omitnan');
beta = (TC_rh_mesh*GS')/(GS*GS');
TC_rh_mesh = TC_rh_mesh - beta*GS;
clearvars GS beta;
%% bandpass filtering
fsTem = 1/0.72;
fLow = 0.01 ;
fHigh = 0.1 ;
isdemean = 1;
% left brain

TC_lh_mesh = bandpa_fMRI(TC_lh_mesh,fsTem,fLow,fHigh,isdemean) ;
TC_lh_mesh_f = vec2map(TC_lh_mesh,cii_mask_left);
% right brain
TC_rh_mesh = bandpa_fMRI(TC_rh_mesh,fsTem,fLow,fHigh,isdemean) ;
TC_rh_mesh_f = vec2map(TC_rh_mesh,cii_mask_right);
clearvars TC_lh_mesh TC_rh_mesh

%zscore
TC_lh_mesh_f = zscore(TC_lh_mesh_f,[],3);
TC_rh_mesh_f = zscore(TC_rh_mesh_f,[],3);

%% save data
data.BrainImg_lh = TC_lh_mesh_f;
data.BrainImg_rh = TC_rh_mesh_f;
data.F_fs = 0.72;
data.F_band = [0.01 1];
data.F_after = 'zscore';%'none' 'demean'
data.R_globalSig = 'yes';

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