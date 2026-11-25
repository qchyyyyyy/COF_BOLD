function preprocessed = meg_bold_preprocess_source(generated, cfg)
%MEG_BOLD_PREPROCESS_SOURCE Apply Project 01 task-fMRI temporal operations.

Y = double(generated.bold_source'); % time x source
Y = bsxfun(@minus,Y,mean(Y,1));
GS = mean(Y,2);
[Yhp,X0] = meg_highpass(Y,cfg.output_tr_seconds,cfg.highpass_seconds);
if isempty(X0), GShp = GS; else, GShp = GS - X0*(X0'*GS); end
denominator = GShp' * GShp;
if denominator <= eps
    error('MEG:WorkflowGlobalSignal', 'High-pass global signal has zero variance.');
end
beta = (GShp' * Yhp) ./ denominator;
residual = Yhp - GShp * beta;

h72 = meg_hrf(cfg.output_tr_seconds,double(cfg.hrf_parameters));
H = fft(h72,size(residual,1));
W = conj(H) ./ (abs(H).^2 + double(cfg.wiener_lambda));
exact = real(ifft(bsxfun(@times,fft(residual,[],1),W),[],1));
exact = block_lowpass(exact,cfg.lowpass_hz,1/cfg.output_tr_seconds);

preprocessed = struct;
preprocessed.exact_project1_source = single(exact');
preprocessed.global_signal_beta = single(beta(:));
preprocessed.highpass_seconds = double(cfg.highpass_seconds);
preprocessed.lowpass_hz = double(cfg.lowpass_hz);
preprocessed.wiener_lambda = double(cfg.wiener_lambda);
preprocessed.output_tr_seconds = double(cfg.output_tr_seconds);
preprocessed.reliable_frame_mask = generated.reliable_frame_mask;
end

function output = block_lowpass(input,cutoff,fs)
if cutoff >= fs/2
    output = input;
    return
end
output = zeros(size(input));
block = 512;
for first = 1:block:size(input,2)
    ix = first:min(first+block-1,size(input,2));
    output(:,ix) = lowpass(input(:,ix),cutoff,fs);
end
end
