function result = meg_fit_streaming_gls(metricDir, V_mask, blockTimes, frameIntervalSeconds, contrasts, options)
%MEG_FIT_STREAMING_GLS Streamed, gap-aware Project-01-style GLS.
%
% The external COF statistics code applies OLS, estimates a per-vertex AR(1)
% coefficient, whitens with a Prais-Winsten transform, and returns a contrast
% cope/varcope/tstat for each spatial vertex.  This implementation keeps the
% same estimand while reading saved DivMap/CurlMap chunks one at a time.  All
% lagged products and whitening are reset at a segment boundary, which is
% required for MEG bad-data gaps.
%
% Inputs
%   metricDir              directory containing divcurl_chunk_*.mat files
%   V_mask                 2-D logical/numeric mask for the 200x500 map grid
%   blockTimes             struct array with start_seconds, end_seconds,
%                          and condition_code fields
%   frameIntervalSeconds   nominal velocity interval length (0.02 s)
%   contrasts              struct array with name, first_code, second_code
%   options                optional struct; include_unlabeled_intervals
%
% Output fields are arrays with dimensions V x contrasts x patterns.

if nargin < 5
    error('MEG:GLSArguments', ...
        'Usage: meg_fit_streaming_gls(metricDir,V_mask,blockTimes,dt,contrasts,options)');
end
if nargin < 6 || isempty(options), options = struct; end
if ~isfield(options, 'include_unlabeled_intervals')
    options.include_unlabeled_intervals = true;
end
if ~isfield(options, 'min_observations'), options.min_observations = 3; end

metricDir = char(metricDir);
if ~isfolder(metricDir)
    error('MEG:GLSMetricsDir', 'Metric directory not found: %s', metricDir);
end
if ndims(V_mask) ~= 2 || ~any(V_mask(:))
    error('MEG:GLSMask', 'V_mask must be a non-empty 2-D mask.');
end
if frameIntervalSeconds <= 0 || ~isfinite(frameIntervalSeconds)
    error('MEG:GLSDt', 'frameIntervalSeconds must be positive and finite.');
end

manifestFile = fullfile(metricDir, 'divcurl_manifest.json');
if ~isfile(manifestFile)
    error('MEG:GLSManifest', 'Missing D/C manifest: %s', manifestFile);
end
manifest = jsondecode(fileread(manifestFile));
if ~isfield(manifest, 'status') || ~strcmpi(char(manifest.status), 'complete')
    error('MEG:GLSManifestStatus', 'D/C manifest is not complete: %s', manifestFile);
end
if ~isfield(manifest, 'chunks') || isempty(manifest.chunks)
    error('MEG:GLSChunks', 'D/C manifest has no chunks.');
end
chunks = manifest.chunks;
[~, order] = sort([chunks.index]);
chunks = chunks(order);
indices = [chunks.index];
if ~isequal(indices, 1:numel(indices))
    error('MEG:GLSChunkIndices', 'D/C chunk indices are not contiguous.');
end

if ~isstruct(blockTimes) || isempty(blockTimes)
    error('MEG:GLSBlocks', 'blockTimes must be a non-empty struct array.');
end
requiredBlockFields = {'start_seconds','end_seconds','condition_code'};
for iField = 1:numel(requiredBlockFields)
    if ~isfield(blockTimes, requiredBlockFields{iField})
        error('MEG:GLSBlocks', 'blockTimes lacks %s.', requiredBlockFields{iField});
    end
end
for iBlock = 1:numel(blockTimes)
    if blockTimes(iBlock).end_seconds <= blockTimes(iBlock).start_seconds || ...
            ~isfinite(blockTimes(iBlock).start_seconds) || ...
            ~isfinite(blockTimes(iBlock).end_seconds) || ...
            ~ismember(blockTimes(iBlock).condition_code, [1 2 3])
        error('MEG:GLSBlocks', 'Invalid block %d.', iBlock);
    end
end

if ~isstruct(contrasts) || isempty(contrasts)
    error('MEG:GLSContrasts', 'contrasts must be a non-empty struct array.');
end
for iContrast = 1:numel(contrasts)
    required = {'name','first_code','second_code'};
    for iField = 1:numel(required)
        if ~isfield(contrasts(iContrast), required{iField})
            error('MEG:GLSContrasts', 'Contrast %d lacks %s.', iContrast, required{iField});
        end
    end
    if contrasts(iContrast).first_code == contrasts(iContrast).second_code || ...
            ~ismember(contrasts(iContrast).first_code, [1 2 3]) || ...
            ~ismember(contrasts(iContrast).second_code, [1 2 3])
        error('MEG:GLSContrasts', 'Invalid condition codes in contrast %d.', iContrast);
    end
end

maskIndex = find(V_mask(:) == 1);
V = numel(maskIndex);
C = numel(contrasts);
Patt = 4;
P = 3;
patternNames = {'Source','Sink','CCVortex','CVortex'};

% Pass 1: OLS sufficient statistics. XTX is contrast-specific but does not
% depend on the spatial map. XTY retains one P x V matrix per contrast/pattern.
XTX_ols = zeros(P, P, C);
XTY_ols = zeros(P, V, C, Patt);
frameCount = zeros(C,1);
conditionCounts = zeros(C,3);
lastFrame = NaN;
lastTime = NaN;
for iChunk = 1:numel(chunks)
    loaded = load_chunk_maps(metricDir, chunks(iChunk), maskIndex, V_mask);
    n = numel(loaded.frame_index);
    ensure_chunk_order(loaded, lastFrame, lastTime, iChunk);
    if n < 1
        error('MEG:GLSEmptyChunk', 'D/C chunk %d contains no intervals.', iChunk);
    end
    starts = segment_starts(loaded.frame_index, lastFrame);
    lastFrame = double(loaded.frame_index(end));
    lastTime = double(loaded.times(end));
    codes = assign_condition_codes(loaded.times, blockTimes, frameIntervalSeconds, ...
        options.include_unlabeled_intervals);
    for iContrast = 1:C
        X = design_for_contrast(codes, contrasts(iContrast));
        XTX_ols(:,:,iContrast) = XTX_ols(:,:,iContrast) + X' * X;
        frameCount(iContrast) = frameCount(iContrast) + n;
        for code = 1:3
            conditionCounts(iContrast, code) = conditionCounts(iContrast, code) + nnz(codes == code);
        end
        for iPattern = 1:Patt
            Y = loaded.patterns(:,:,iPattern)'; % n x V
            XTY_ols(:,:,iContrast,iPattern) = ...
                XTY_ols(:,:,iContrast,iPattern) + X' * Y;
        end
    end
    clear loaded starts codes
end

if any(frameCount < options.min_observations)
    error('MEG:GLSObservations', 'Insufficient observations in one or more contrasts.');
end
for iContrast = 1:C
    if rank(XTX_ols(:,:,iContrast)) < P || ...
            any(conditionCounts(iContrast, contrasts(iContrast).first_code) == 0) || ...
            any(conditionCounts(iContrast, contrasts(iContrast).second_code) == 0)
        error('MEG:GLSDesignRank', 'Contrast %s has a rank-deficient or empty design.', ...
            char(contrasts(iContrast).name));
    end
end

betaOLS = zeros(P, V, C, Patt);
for iContrast = 1:C
    for iPattern = 1:Patt
        betaOLS(:,:,iContrast,iPattern) = XTX_ols(:,:,iContrast) \ ...
            XTY_ols(:,:,iContrast,iPattern);
    end
end

% Pass 2: estimate a vertex/pattern/contrast-specific AR(1) coefficient from
% residual pairs that remain within a continuous valid segment.
rhoNumerator = zeros(V, C, Patt);
rhoDenominator = zeros(V, C, Patt);
previousResidual = zeros(V, C, Patt);
havePreviousResidual = false;
lastFrame = NaN;
lastTime = NaN;
for iChunk = 1:numel(chunks)
    loaded = load_chunk_maps(metricDir, chunks(iChunk), maskIndex, V_mask);
    n = numel(loaded.frame_index);
    ensure_chunk_order(loaded, lastFrame, lastTime, iChunk);
    starts = segment_starts(loaded.frame_index, lastFrame);
    lastFrame = double(loaded.frame_index(end));
    lastTime = double(loaded.times(end));
    codes = assign_condition_codes(loaded.times, blockTimes, frameIntervalSeconds, ...
        options.include_unlabeled_intervals);
    for iContrast = 1:C
        X = design_for_contrast(codes, contrasts(iContrast));
        for iPattern = 1:Patt
            Y = loaded.patterns(:,:,iPattern)'; % n x V
            E = Y - X * betaOLS(:,:,iContrast,iPattern);
            [rhoNumerator(:,iContrast,iPattern), rhoDenominator(:,iContrast,iPattern)] = ...
                add_residual_pairs(rhoNumerator(:,iContrast,iPattern), ...
                rhoDenominator(:,iContrast,iPattern), E, starts, ...
                previousResidual(:,iContrast,iPattern), havePreviousResidual);
            previousResidual(:,iContrast,iPattern) = E(end,:)';
        end
    end
    havePreviousResidual = true;
    clear loaded starts codes
end
rhoHat = rhoNumerator ./ (rhoDenominator + eps);
rhoHat(~isfinite(rhoHat)) = 0;
rhoHat = max(min(rhoHat, 0.99), -0.99);

% Pass 3: accumulate the per-vertex whitened normal equations.  Expanding
% (x_t-r*x_{t-1}) terms avoids a vertex loop for every frame while retaining
% the exact per-vertex AR(1) whitening used by the external implementation.
Gxx = zeros(P, P, V, C, Patt);
Gxy = zeros(P, V, C, Patt);
Gy2 = zeros(V, C, Patt);
nWhitened = zeros(C,1);
lastFrame = NaN;
lastTime = NaN;
previousY = zeros(V, Patt);
previousX = zeros(1, P, C);
havePreviousY = false;
for iChunk = 1:numel(chunks)
    loaded = load_chunk_maps(metricDir, chunks(iChunk), maskIndex, V_mask);
    n = numel(loaded.frame_index);
    ensure_chunk_order(loaded, lastFrame, lastTime, iChunk);
    starts = segment_starts(loaded.frame_index, lastFrame);
    lastFrame = double(loaded.frame_index(end));
    lastTime = double(loaded.times(end));
    codes = assign_condition_codes(loaded.times, blockTimes, frameIntervalSeconds, ...
        options.include_unlabeled_intervals);
    for iContrast = 1:C
        X = design_for_contrast(codes, contrasts(iContrast));
        for iPattern = 1:Patt
            Y = loaded.patterns(:,:,iPattern)'; % n x V
            rho = rhoHat(:,iContrast,iPattern);
            % MATLAB does not permit an indexed expression as a function
            % output.  Use local temporaries, then write the updated slices
            % back to the accumulator arrays.
            gxxTmp = Gxx(:,:,:,iContrast,iPattern);
            gxyTmp = Gxy(:,:,iContrast,iPattern);
            gy2Tmp = Gy2(:,iContrast,iPattern);
            [gxxTmp, gxyTmp, gy2Tmp] = add_whitened_chunk( ...
                gxxTmp, gxyTmp, gy2Tmp, X, Y, starts, rho, ...
                previousY(:,iPattern), previousX(1,:,iContrast), havePreviousY);
            Gxx(:,:,:,iContrast,iPattern) = gxxTmp;
            Gxy(:,:,iContrast,iPattern) = gxyTmp;
            Gy2(:,iContrast,iPattern) = gy2Tmp;
        end
        previousX(1,:,iContrast) = X(end,:);
        nWhitened(iContrast) = nWhitened(iContrast) + n;
    end
    % The lagged response state is shared across contrasts, but must remain
    % the previous chunk's final observation until every contrast has used it.
    % Updating it inside the contrast loop would make contrast 2/3 whiten the
    % first row against the current (rather than previous) chunk endpoint.
    for iPattern = 1:Patt
        previousY(:,iPattern) = loaded.patterns(:,end,iPattern);
    end
    havePreviousY = true;
    clear loaded starts codes
end

cope = nan(V, C, Patt);
varcope = nan(V, C, Patt);
tstat = nan(V, C, Patt);
betaGLS = zeros(P, V, C, Patt);
c = [0; 1; -1];
rankDesign = P;
for iContrast = 1:C
    for iPattern = 1:Patt
        for v = 1:V
            thisG = Gxx(:,:,v,iContrast,iPattern);
            thisG = (thisG + thisG') ./ 2;
            thisB = thisG \ Gxy(:,v,iContrast,iPattern);
            betaGLS(:,v,iContrast,iPattern) = thisB;
            sse = Gy2(v,iContrast,iPattern) - 2 .* (thisB' * Gxy(:,v,iContrast,iPattern)) + ...
                thisB' * thisG * thisB;
            dof = max(nWhitened(iContrast) - rankDesign, 1);
            sig2 = max(0, sse ./ dof);
            cCov = (c' / thisG) * c;
            cope(v,iContrast,iPattern) = c' * thisB;
            varcope(v,iContrast,iPattern) = sig2 .* cCov;
            tstat(v,iContrast,iPattern) = cope(v,iContrast,iPattern) ./ ...
                sqrt(varcope(v,iContrast,iPattern) + eps);
        end
    end
end

result = struct;
result.schema_version = 1;
result.solver = 'fit_gls_ar1';
result.segment_aware = true;
result.metric_manifest = manifestFile;
result.V_mask = V_mask;
result.mask_index = maskIndex;
result.pattern_names = patternNames;
result.contrasts = contrasts;
result.contrast_vector = c;
result.cope = cope;
result.varcope = varcope;
result.tstat = tstat;
result.beta_gls = betaGLS;
result.beta_ols = betaOLS;
result.rho_hat = rhoHat;
result.frame_count = frameCount;
result.condition_counts = conditionCounts;
result.whitened_frame_count = nWhitened;
result.options = options;
end

function loaded = load_chunk_maps(metricDir, chunk, maskIndex, V_mask)
filename = fullfile(metricDir, char(chunk.filename));
if ~isfile(filename)
    error('MEG:GLSChunkMissing', 'Missing D/C map chunk: %s', filename);
end
loadedRaw = load(filename, 'DivMap', 'CurlMap', ...
    'interval_frame_index', 'interval_time_seconds');
required = {'DivMap','CurlMap','interval_frame_index','interval_time_seconds'};
for i = 1:numel(required)
    if ~isfield(loadedRaw, required{i})
        error('MEG:GLSChunkFields', 'Chunk %s lacks %s.', filename, required{i});
    end
end
n = numel(loadedRaw.interval_frame_index);
if n < 1
    error('MEG:GLSEmptyChunk', 'D/C chunk contains no intervals: %s', filename);
end
frameIndex = double(loadedRaw.interval_frame_index(:)');
times = double(loadedRaw.interval_time_seconds(:)');
if numel(times) ~= n || any(~isfinite(frameIndex)) || any(~isfinite(times)) || ...
        (n > 1 && (any(diff(frameIndex) <= 0) || any(diff(times) <= 0)))
    error('MEG:GLSChunkOrder', ...
        'Frame indices/times must be finite and strictly increasing in %s.', filename);
end
if size(loadedRaw.DivMap,1) ~= size(V_mask,1) || ...
        size(loadedRaw.DivMap,2) ~= size(V_mask,2) || ...
        size(loadedRaw.DivMap,3) ~= n || ...
        size(loadedRaw.CurlMap,1) ~= size(V_mask,1) || ...
        size(loadedRaw.CurlMap,2) ~= size(V_mask,2) || ...
        size(loadedRaw.CurlMap,3) ~= n || ...
        ~isequal(size(loadedRaw.DivMap), size(loadedRaw.CurlMap))
    error('MEG:GLSChunkShape', 'Map/frame shape mismatch in %s.', filename);
end
div = reshape(double(loadedRaw.DivMap), [], n);
curl = reshape(double(loadedRaw.CurlMap), [], n);
div = div(maskIndex,:);
curl = curl(maskIndex,:);
if any(~isfinite(div(:))) || any(~isfinite(curl(:)))
    error('MEG:GLSChunkFinite', 'Nonfinite map value inside V_mask in %s.', filename);
end
% Occupancy maps are nonnegative and retain the Project-01 Source/Sink and
% CCVortex/CVortex polarity convention.
patterns = cat(3, div > 0, div < 0, curl > 0, curl < 0);
loaded = struct('patterns', double(patterns), ...
    'frame_index', frameIndex, 'times', times);
end

function starts = segment_starts(frameIndex, previousFrame)
frameIndex = double(frameIndex(:)');
starts = false(1, numel(frameIndex));
if isempty(frameIndex), return; end
starts(1) = isnan(previousFrame) || frameIndex(1) ~= previousFrame + 1;
if numel(frameIndex) > 1
    starts(2:end) = diff(frameIndex) ~= 1;
end
end

function ensure_chunk_order(loaded, previousFrame, previousTime, chunkIndex)
if isempty(loaded.frame_index) || isempty(loaded.times)
    error('MEG:GLSEmptyChunk', 'D/C chunk %d contains no intervals.', chunkIndex);
end
if ~isnan(previousFrame) && (loaded.frame_index(1) <= previousFrame || ...
        loaded.times(1) <= previousTime)
    error('MEG:GLSChunkOrder', ...
        'D/C chunks are not strictly ordered at chunk %d.', chunkIndex);
end
end

function codes = assign_condition_codes(times, blockTimes, dt, includeUnlabeled)
times = double(times(:)');
codes = zeros(size(times));
tol = max(1e-9, dt * 1e-4);
for iBlock = 1:numel(blockTimes)
    in = times >= blockTimes(iBlock).start_seconds - tol & ...
        times + dt <= blockTimes(iBlock).end_seconds + tol;
    if any(in & codes ~= 0)
        error('MEG:GLSBlockOverlap', 'An interval is assigned to multiple blocks.');
    end
    codes(in) = blockTimes(iBlock).condition_code;
end
if ~includeUnlabeled
    codes(codes == 0) = NaN;
end
end

function X = design_for_contrast(codes, contrast)
codes = double(codes(:));
if any(~isfinite(codes))
    keep = isfinite(codes);
    codes = codes(keep);
else
    keep = true(size(codes));
end
% This helper is called with finite code vectors in the current stream. The
% keep branch is retained so a future caller can explicitly exclude padding.
X = [ones(numel(codes),1), double(codes == contrast.first_code), ...
    double(codes == contrast.second_code)];
if ~all(keep)
    error('MEG:GLSNaNDesign', 'NaN condition labels require synchronized data filtering.');
end
end

function [numOut, denOut] = add_residual_pairs(numIn, denIn, E, starts, previousE, havePrevious)
numOut = numIn;
denOut = denIn;
if isempty(E), return; end
% Cross-chunk pair, only when the first frame continues the previous segment.
if havePrevious && ~starts(1)
    ePrev = previousE(:);
    eCurr = E(1,:)';
    numOut = numOut + ePrev .* eCurr;
    denOut = denOut + ePrev .* ePrev;
end
if size(E,1) > 1
    pairKeep = ~starts(2:end);
    if any(pairKeep)
        prev = E(1:end-1,:); curr = E(2:end,:);
        prev = prev(pairKeep,:); curr = curr(pairKeep,:);
        numOut = numOut + sum(prev .* curr, 1)';
        denOut = denOut + sum(prev .* prev, 1)';
    end
end
end

function [GxxOut, GxyOut, Gy2Out] = add_whitened_chunk(GxxIn, GxyIn, Gy2In, ...
        X, Y, starts, rho, previousY, previousX, havePrevious)
% Add start rows and within-segment continuation rows using the expanded
% Prais-Winsten equations. Y is n x V; rho is V x 1.
P = size(X,2);
V = size(Y,2);
GxxOut = GxxIn;
GxyOut = GxyIn;
Gy2Out = Gy2In;
startRows = find(starts);
if ~isempty(startRows)
    Xs = X(startRows,:);
    A = Xs' * Xs;
    Sxy = Xs' * Y(startRows,:);
    Sy2 = sum(Y(startRows,:).^2, 1)';
    oneMinus = 1 - rho.^2;
    for a = 1:P
        for b = 1:P
            GxxOut(a,b,:) = GxxOut(a,b,:) + reshape(A(a,b) .* oneMinus, 1, 1, V);
        end
        % Both the whitened design and response carry the Prais--Winsten
        % sqrt(1-rho^2) factor, so their cross-product receives 1-rho^2.
        GxyOut(a,:) = GxyOut(a,:) + Sxy(a,:) .* oneMinus';
    end
    Gy2Out = Gy2Out + Sy2 .* oneMinus;
end

% Within-chunk continuation rows.
pairRows = find(~starts(2:end)) + 1;
if ~isempty(pairRows)
    Xcurr = X(pairRows,:);
    Xprev = X(pairRows-1,:);
    Ycurr = Y(pairRows,:);
    Yprev = Y(pairRows-1,:);
    [GxxOut, GxyOut, Gy2Out] = add_continuation(...
        Xcurr, Xprev, Ycurr, Yprev, rho, GxxOut, GxyOut, Gy2Out);
end

% The first row of a chunk can continue the previous chunk/segment.
if havePrevious && ~starts(1)
    Xcurr = X(1,:); Xprev = previousX;
    Ycurr = Y(1,:); Yprev = previousY(:)';
    [GxxOut, GxyOut, Gy2Out] = add_continuation(...
        Xcurr, Xprev, Ycurr, Yprev, rho, GxxOut, GxyOut, Gy2Out);
end
end

function [GxxOut, GxyOut, Gy2Out] = add_continuation(Xcurr, Xprev, Ycurr, Yprev, ...
        rho, GxxOut, GxyOut, Gy2Out)
% Add continuation observations. This helper accepts either many rows or one.
if size(Xcurr,1) == 1
    A = Xcurr' * Xcurr;
    B = Xcurr' * Xprev;
    C = Xprev' * Xprev;
    Dcc = Xcurr' .* Ycurr;
    Dcp = Xcurr' .* Yprev;
    Dpc = Xprev' .* Ycurr;
    Dpp = Xprev' .* Yprev;
    y2 = Ycurr(:).^2;
    crossY = Yprev(:) .* Ycurr(:);
    prevY2 = Yprev(:).^2;
else
    A = Xcurr' * Xcurr;
    B = Xcurr' * Xprev;
    C = Xprev' * Xprev;
    Dcc = Xcurr' * Ycurr;
    Dcp = Xcurr' * Yprev;
    Dpc = Xprev' * Ycurr;
    Dpp = Xprev' * Yprev;
    y2 = sum(Ycurr.^2,1)';
    crossY = sum(Yprev .* Ycurr,1)';
    prevY2 = sum(Yprev.^2,1)';
end
rhoRow = rho(:)';
rho2Row = (rho.^2)';
for a = 1:size(Xcurr,2)
    for b = 1:size(Xcurr,2)
        GxxOut(a,b,:) = GxxOut(a,b,:) + reshape( ...
            A(a,b) - rhoRow .* (B(a,b) + B(b,a)) + rho2Row .* C(a,b), 1, 1, []);
    end
    % Full Prais--Winsten cross-product:
    % (X_t-rho X_{t-1})' (Y_t-rho Y_{t-1}).  Keeping all four terms is
    % necessary for numerical equivalence to the external fit_gls_ar1.m.
    GxyOut(a,:) = GxyOut(a,:) + Dcc(a,:) - Dcp(a,:) .* rhoRow - ...
        Dpc(a,:) .* rhoRow + Dpp(a,:) .* rho2Row;
end
Gy2Out = Gy2Out + y2 - 2 .* rho .* crossY + (rho.^2) .* prevY2;
end
