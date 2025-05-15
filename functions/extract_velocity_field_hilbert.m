function [Ux,Uy] = extract_velocity_field_hilbert(BrainImg,mask,rho,tau,sigma)
% Input:
%   path: data storage address
%   brainimg_filename: filename of preprocessed fMRI data. The pipeline of
%   preprocessing incluedes...
%       1. global signal regression
%       2. interpolation 
%       3. bandpass filter (0.01-0.1Hz)
%   hemisphrere: hemisphrere of cerebral cortex ('lh' or 'rh')
%   mask: a mask for image
%   rho: the penalty parameter for the smooth regularization
%   tau: the penalty parameter for the uniqueness of the optimal solution
%   sigma: the parameter for detecting VVFs at certain spatial scale
% Output:
%   Ux: 
%   Uy: 
%
%   by Chengyuan Qian @Fudan November 14, 2023

% ---------- Convert data to analytic data via hilbert transformation ----------
if isstruct(BrainImg)
    load([BrainImg.path BrainImg.fileName],'data');
    hemisphere = BrainImg.hemisphere;
    savePath = BrainImg.path;
    switch hemisphere
        case 'lh'
            BrainImg = data.BrainImg_lh;
        case 'rh'
            BrainImg = data.BrainImg_rh;
    end
    data = rmfield(data,'BrainImg_lh');
    data = rmfield(data,'BrainImg_rh');
else
    savePath='';
end

%% hilbert transformation
BrainImg = map2vec(BrainImg,mask);
if isreal(BrainImg)
    BrainImg = hilbert(BrainImg').';
end
%%

% ---------- Create an 1*2 cell that contains the real and imaginary parts of 
%the analytic BOLD-fMRI signal ((M+1)*(N+1)*T tensor y-, x-, t-direction) ----------  

complex_Img{1}=vec2map(real(BrainImg),mask);
complex_Img{2} =vec2map(imag(BrainImg),mask);
clearvars BrainImg
% ---------- Calculate velocity vector field of analytic signal ----------
para_model.tau = tau;
para_model.rho = rho;
para_model.sigma = sigma;

[Ux,Uy] = velocity_field_multidim_ell2(complex_Img,mask,para_model);

% ---------- Save results ----------
if exist(savePath,'dir')
    VelField.Ux = Ux;
    VelField.Uy = Uy;
    VelField.para_model = para_model;
    VelField.data = data;

    save([savePath,'VelField_',hemisphere,'_hilbert.mat'],'VelField','-v7.3');
end