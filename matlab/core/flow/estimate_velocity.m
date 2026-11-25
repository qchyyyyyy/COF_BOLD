function [Ux,Uy] = estimate_velocity(BrainImg,mask,method,rho,tau,sigma,boundaryExtensionRings)
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
%   boundaryExtensionRings: optional harmonic-fill rings for unsmoothed COF
% Output:
%   Ux:
%   Uy:
%
%   by Chengyuan Qian @Fudan November 14, 2023

% ---------- Convert data to analytic data via hilbert transformation ----------

assert(ismember(string(method),["COF","OF","PG"]),'cof:Method','Unknown velocity estimator.');
validateattributes(rho,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(tau,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(sigma,{'numeric'},{'scalar','nonnegative','finite'});
if nargin < 7
    boundaryExtensionRings = [];
end

if isstruct(BrainImg)
    hemisphere = BrainImg.hemisphere;
    assert(ismember(string(hemisphere),["lh","rh"]),'cof:Hemisphere','Expected lh or rh.');
    savePath = BrainImg.path;

    % if exist([savePath,'VelField_',hemisphere,'_' method '.mat'],'file')
    %     disp(['VelField_',hemisphere,'_' method '.mat already exist']);
    %     return
    % end

    inputFile = fullfile(BrainImg.path, BrainImg.fileName);
    if exist(inputFile,'file')
        load(inputFile,'data');
    else
        error('estimate_velocity:MissingInput', ...
            'Preprocessed input does not exist: %s', inputFile);
    end

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

switch method
    case 'COF'
        BrainImg = grid_to_vector(BrainImg,mask);
        if isreal(BrainImg)
            BrainImg = hilbert(BrainImg').';
        end

        % ---------- Create an 1*2 cell that contains the real and imaginary parts of
        %the analytic BOLD-fMRI signal ((M+1)*(N+1)*T tensor y-, x-, t-direction) ----------

        complexImg{1}=vector_to_grid(real(BrainImg),mask);
        complexImg{2} =vector_to_grid(imag(BrainImg),mask);
        clearvars BrainImg
        % ---------- Calculate velocity vector field of analytic signal ----------
        para_model.tau = tau;
        para_model.rho = rho;
        para_model.sigma = sigma;
        if ~isempty(boundaryExtensionRings)
            para_model.boundaryExtensionRings = boundaryExtensionRings;
        end
        [Ux,Uy] = cof_solve(complexImg,mask,para_model);

    case 'PG'
        filtSigma = sigma;
        filtWidth = ceil(3*filtSigma);
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        BrainImg(repmat(mask,1,1,size(BrainImg,3))==0)=nan;
        for t = 1:size(BrainImg,3)
            BrainImg(:,:,t) = nanconv(BrainImg(:,:,t),imageFilter);
        end
        mask0 = double(~isnan(BrainImg(:,:,1)));
        BrainImg = grid_to_vector(BrainImg,mask0);
        if isreal(BrainImg)
            phaseImg = angle(hilbert(BrainImg').');
        else
            phaseImg = angle(BrainImg);
        end
        phaseImg = vector_to_grid(phaseImg,mask0);

        [Ux,Uy] = phase_gradient(phaseImg);
    case 'OF'
        para_model.tau = tau;
        para_model.rho = rho;
        para_model.sigma = sigma;
        [Ux,Uy] = of_solve(BrainImg,mask,para_model);
end

% ---------- Save results ----------
if exist(savePath,'dir')
    VelField.Ux = Ux;
    VelField.Uy = Uy;
    % Preserve the historical output contract: all estimators except PG
    % carry their parameter struct in the saved velocity-field file.
    if ismember(string(method),["COF","OF"]), VelField.para_model = para_model; end
    VelField.data = data;

    save(fullfile(savePath, ['VelField_',hemisphere,'_' method '.mat']), ...
        'VelField','-v7.3');
end
