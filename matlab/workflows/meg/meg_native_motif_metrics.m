function manifest = meg_native_motif_metrics(cofDir, outputDir, chunkRange, options)
%MEG_NATIVE_MOTIF_METRICS Compute D/C maps and global PQ/PS per frame.
%
% The divergence/curl maps are calculated chunk-wise with meg_velocity_derivatives.m.
% Their z-score and +/- n_std threshold are estimated over all COF chunks,
% not independently per chunk. chunkRange is useful for parallel execution;
% every worker still performs the same full-session first pass for statistics.

if nargin < 3 || isempty(chunkRange), chunkRange = [1 Inf]; end
if nargin < 4 || isempty(options), options = struct; end
if ~(ischar(cofDir) || isstring(cofDir)) || ~(ischar(outputDir) || isstring(outputDir))
    error('MEG:MetricsArguments', 'cofDir and outputDir must be paths.');
end
cofDir = char(cofDir);
outputDir = char(outputDir);
if numel(chunkRange) ~= 2 || chunkRange(1) < 1 || chunkRange(2) < chunkRange(1)
    error('MEG:MetricsChunkRange', 'chunkRange must be [first last].');
end
if ~isfield(options,'filtSigma'), options.filtSigma = 3; end
if ~isfield(options,'doZscore'), options.doZscore = true; end
if ~isfield(options,'doBinarize'), options.doBinarize = true; end
if ~isfield(options,'nStd'), options.nStd = 2; end
if ~isfield(options,'sizeThreshold'), options.sizeThreshold = 10; end
if ~isfield(options,'saveMaps'), options.saveMaps = true; end

cofManifest = jsondecode(fileread(fullfile(cofDir,'cof_manifest.json')));
chunks = cofManifest.velocity_chunks;
nChunks = numel(chunks);
firstChunk = min(max(1,chunkRange(1)), nChunks);
lastChunk = min(chunkRange(2), nChunks);
loadedPreprocessing = load(fullfile(cofDir,'cof_preprocessing.mat'), 'preprocessing');
mask = loadedPreprocessing.preprocessing.mask;
V_mask = conv2(mask, ones(2,2)./4, 'valid');
V_mask(V_mask ~= 1) = 0;

if isfolder(outputDir)
    existing = dir(outputDir);
    existing = existing(~ismember({existing.name},{'.','..'}));
    if ~isempty(existing)
        error('MEG:MetricsExistingOutput', 'Output directory is not empty: %s', outputDir);
    end
else
    mkdir(outputDir);
end
write_json(fullfile(outputDir,'divcurl_progress.json'), ...
    struct('status','statistics_pass','chunks_total',nChunks, ...
    'chunks_completed',0,'chunks_written',0));

% Pass 1: session-global raw divergence/curl moments.
sumDiv = 0; sumSqDiv = 0; countDiv = 0;
sumCurl = 0; sumSqCurl = 0; countCurl = 0;
for iChunk = 1:nChunks
    loaded = load(fullfile(cofDir,chunks(iChunk).filename), 'Ux','Uy');
    [divMap,curlMap] = meg_velocity_derivatives(loaded.Ux, loaded.Uy, V_mask, ...
        options.filtSigma, 0, 0, options.nStd);
    finiteDiv = isfinite(divMap);
    finiteCurl = isfinite(curlMap);
    divValues = double(divMap(finiteDiv));
    curlValues = double(curlMap(finiteCurl));
    sumDiv = sumDiv + sum(divValues);
    sumSqDiv = sumSqDiv + sum(divValues.^2);
    countDiv = countDiv + numel(divValues);
    sumCurl = sumCurl + sum(curlValues);
    sumSqCurl = sumSqCurl + sum(curlValues.^2);
    countCurl = countCurl + numel(curlValues);
    if mod(iChunk,10) == 0 || iChunk == nChunks
        write_json(fullfile(outputDir,'divcurl_progress.json'), ...
            struct('status','statistics_pass','chunks_total',nChunks, ...
            'chunks_completed',iChunk,'chunks_written',0));
    end
    clear loaded divMap curlMap finiteDiv finiteCurl divValues curlValues
end
divMean = sumDiv / countDiv;
curlMean = sumCurl / countCurl;
if countDiv <= 1 || countCurl <= 1
    error('MEG:MetricsVariance', 'Insufficient finite divergence/curl values.');
end
divStd = sqrt(max(0, (sumSqDiv - sumDiv.^2./countDiv) ./ (countDiv-1)));
curlStd = sqrt(max(0, (sumSqCurl - sumCurl.^2./countCurl) ./ (countCurl-1)));
if options.doZscore && (divStd <= eps || curlStd <= eps)
    error('MEG:MetricsVariance', 'Divergence or curl has zero session variance.');
end

stats = struct('filtSigma',options.filtSigma,'do_zscore',options.doZscore, ...
    'do_binarize',options.doBinarize,'n_std',options.nStd, ...
    'size_threshold',options.sizeThreshold,'div_mean',divMean, ...
    'div_std',divStd,'curl_mean',curlMean,'curl_std',curlStd, ...
    'finite_div_values',countDiv,'finite_curl_values',countCurl);
save(fullfile(outputDir,'divcurl_statistics.mat'),'stats','V_mask');

chunkInfo = struct('index',{},'filename',{},'interval_start',{}, ...
    'interval_end',{},'frames',{},'bytes',{});
for iChunk = firstChunk:lastChunk
    loaded = load(fullfile(cofDir,chunks(iChunk).filename), 'Ux','Uy', ...
        'interval_frame_index','interval_time_seconds');
    [divMap,curlMap] = meg_velocity_derivatives(loaded.Ux, loaded.Uy, V_mask, ...
        options.filtSigma, 0, 0, options.nStd);
    interval_frame_index = loaded.interval_frame_index;
    interval_time_seconds = loaded.interval_time_seconds;
    if options.doZscore
        divMap = (divMap - divMean) ./ divStd;
        curlMap = (curlMap - curlMean) ./ curlStd;
    end
    if options.doBinarize
        if options.doZscore
            divMap(divMap < options.nStd & divMap > -options.nStd) = 0;
            curlMap(curlMap < options.nStd & curlMap > -options.nStd) = 0;
        else
            divMap(divMap < divMean+options.nStd*divStd & divMap > divMean-options.nStd*divStd) = 0;
            curlMap(curlMap < curlMean+options.nStd*curlStd & curlMap > curlMean-options.nStd*curlStd) = 0;
        end
        divMap(divMap > 0) = 1; divMap(divMap < 0) = -1;
        curlMap(curlMap > 0) = 1; curlMap(curlMap < 0) = -1;
    end
    nIntervals = size(divMap,3);
    NumSource = zeros(1,nIntervals);
    NumSink = zeros(1,nIntervals);
    NumCCVortex = zeros(1,nIntervals);
    NumCVortex = zeros(1,nIntervals);
    SizeSource = zeros(1,nIntervals);
    SizeSink = zeros(1,nIntervals);
    SizeCCVortex = zeros(1,nIntervals);
    SizeCVortex = zeros(1,nIntervals);
    for iTime = 1:nIntervals
        [NumSource(iTime),NumSink(iTime)] = meg_count_motifs(divMap(:,:,iTime), options.sizeThreshold);
        [NumCCVortex(iTime),NumCVortex(iTime)] = meg_count_motifs(curlMap(:,:,iTime), options.sizeThreshold);
        [SizeSource(iTime),SizeSink(iTime)] = meg_measure_motif_sizes(divMap(:,:,iTime));
        [SizeCCVortex(iTime),SizeCVortex(iTime)] = meg_measure_motif_sizes(curlMap(:,:,iTime));
    end
    PQ = NumSource + NumSink;
    PS = SizeSource + SizeSink;
    PQ_curl = NumCCVortex + NumCVortex;
    PS_curl = SizeCCVortex + SizeCVortex;
    filename = sprintf('divcurl_chunk_%04d.mat',iChunk);
    if options.saveMaps
        DivMap = single(divMap); %#ok<NASGU>
        CurlMap = single(curlMap); %#ok<NASGU>
        save(fullfile(outputDir,filename),'DivMap','CurlMap','NumSource', ...
            'NumSink','NumCCVortex','NumCVortex','SizeSource','SizeSink', ...
            'SizeCCVortex','SizeCVortex','PQ','PS','PQ_curl','PS_curl', ...
            'interval_frame_index','interval_time_seconds','stats','-v7.3');
    else
        save(fullfile(outputDir,filename),'NumSource','NumSink','NumCCVortex', ...
            'NumCVortex','SizeSource','SizeSink','SizeCCVortex','SizeCVortex', ...
            'PQ','PS','PQ_curl','PS_curl','interval_frame_index', ...
            'interval_time_seconds','stats','-v7.3');
    end
    info = dir(fullfile(outputDir,filename));
    chunkInfo(end+1).index = iChunk; %#ok<AGROW>
    chunkInfo(end).filename = filename;
    chunkInfo(end).interval_start = chunks(iChunk).interval_start;
    chunkInfo(end).interval_end = chunks(iChunk).interval_end;
    chunkInfo(end).frames = nIntervals;
    chunkInfo(end).bytes = info.bytes;
    progress = struct('status','running','chunk_start',firstChunk, ...
        'chunk_end',lastChunk,'chunks_written',numel(chunkInfo), ...
        'last_chunk',iChunk,'intervals_written',sum([chunkInfo.frames]));
    write_json(fullfile(outputDir,'divcurl_progress.json'),progress);
end

manifest = struct;
manifest.schema_version = 1;
manifest.status = 'complete';
manifest.input_cof_manifest = fullfile(cofDir,'cof_manifest.json');
manifest.method = 'meg_velocity_derivatives.m';
manifest.pattern_functions = {'meg_count_motifs.m','meg_measure_motif_sizes.m'};
manifest.chunk_range = [firstChunk,lastChunk];
manifest.chunk_count = numel(chunkInfo);
manifest.parameters = stats;
manifest.velocity_parameters = cofManifest.parameters;
manifest.metrics_convention = 'PQ=NumSource+NumSink; PS=SizeSource+SizeSink; curl variants retained separately';
manifest.chunks = chunkInfo;
write_json(fullfile(outputDir,'divcurl_manifest.json'),manifest);
write_json(fullfile(outputDir,'divcurl_progress.json'), ...
    struct('status','complete','chunk_start',firstChunk,'chunk_end',lastChunk, ...
    'chunks_written',numel(chunkInfo),'last_chunk',lastChunk, ...
    'intervals_written',sum([chunkInfo.frames])));
end

function write_json(filename,value)
fid = fopen(filename,'w');
if fid < 0, error('MEG:MetricsManifest', 'Cannot open %s.',filename); end
cleanup = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(value),'char');
clear cleanup
end
