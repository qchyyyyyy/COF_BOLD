function generated = meg_bold_generate(power, validMask, cfg)
%MEG_BOLD_GENERATE Convert source power to a 720-ms HRF observation.

power = single(power);
validMask = logical(validMask(:)');
[nSources,nFrames] = size(power);
if numel(validMask) ~= nFrames || nnz(validMask) < 2
    error('MEG:WorkflowPowerShape', 'Power and valid-frame mask are inconsistent.');
end
nativeDt = double(cfg.native_dt_seconds);
outputTr = double(cfg.output_tr_seconds);
step = round(outputTr/nativeDt);
h = meg_hrf(nativeDt, double(cfg.hrf_parameters));
nFft = 2^nextpow2(nFrames + numel(h) - 1);
H = fft(double(h(:)'),nFft,2);
outputIndex = 1:step:nFrames;
nOutput = numel(outputIndex);
validIndex = find(validMask);
firstValid = validIndex(1); lastValid = validIndex(end);

sourceMean = mean(double(power(:,validMask)),2);
if any(~isfinite(sourceMean) | sourceMean <= eps)
    error('MEG:WorkflowSourceMean', 'Every source must have positive finite mean power.');
end
bold = zeros(nSources,nOutput,'single');
driveDownsampled = zeros(nSources,nOutput,'single');
blockSize = double(cfg.source_block_size);
for first = 1:blockSize:nSources
    ix = first:min(first+blockSize-1,nSources);
    x = double(power(ix,:)) ./ sourceMean(ix) - 1;
    % Zero relative power is the registered history/edge assumption.  Only
    % internal invalid samples are linearly interpolated.
    filled = zeros(size(x));
    filled(:,validMask) = x(:,validMask);
    internal = firstValid:lastValid;
    if any(~validMask(internal))
        filled(:,internal) = interp1(double(validIndex),x(:,validMask)', ...
            double(internal),'linear')';
    end
    spectrum = fft(filled,nFft,2);
    convolved = real(ifft(bsxfun(@times,spectrum,H),[],2));
    convolved = convolved(:,1:nFrames);
    bold(ix,:) = single(convolved(:,outputIndex));
    driveDownsampled(ix,:) = single(filled(:,outputIndex));
end

gapInfluence = filter(abs(double(h(:)')),1,double(~validMask));
absoluteMass = sum(abs(double(h)));
gapInfluence = gapInfluence ./ absoluteMass;
edgeInfluence = zeros(1,nFrames);
tailMass = flip(cumsum(flip(abs(double(h(:))))));
for t = 1:min(nFrames,numel(h)-1)
    edgeInfluence(t) = tailMass(t+1) ./ absoluteMass;
end
imputedInfluence = min(1,gapInfluence + edgeInfluence);
outputInfluence = imputedInfluence(outputIndex);
reliable = outputInfluence <= double(cfg.maximum_imputed_hrf_mass);

generated = struct;
generated.schema_version = 1;
generated.estimand = 'canonical-HRF observation of relative fixed-orientation LCMV power';
generated.bold_source = bold;
generated.drive_downsampled = driveDownsampled;
generated.source_mean_power = single(sourceMean);
generated.output_time_seconds = (outputIndex-1).*nativeDt;
generated.output_native_frame_index = outputIndex;
generated.imputed_hrf_mass_fraction = single(outputInfluence);
generated.reliable_frame_mask = reliable;
generated.native_valid_frame_count = nnz(validMask);
generated.native_invalid_frame_count = nnz(~validMask);
generated.internal_gap_policy = 'per-source linear interpolation';
generated.edge_history = 'zero relative power';
generated.native_dt_seconds = nativeDt;
generated.output_tr_seconds = outputTr;
generated.hrf = single(h);
generated.hrf_parameters = double(cfg.hrf_parameters(:)');
end
