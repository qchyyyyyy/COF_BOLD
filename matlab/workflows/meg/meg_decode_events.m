function decoded = meg_decode_events(rawFile, megconnectomeRoot, preStimTime, postStimTime)
%MEG_DECODE_EVENTS Decode official HCP Working-Memory events.
%
% The official trial function is called in three modes.  Block and fixation
% calls provide exact raw-sample boundaries; the trial call provides the
% official image-locked task windows and Fix15sec pseudo-trials.  No labels are
% inferred from the already processed PQ/PS derivatives.

if nargin < 2 || isempty(megconnectomeRoot)
    error('MEG:EventDecoderArguments', 'rawFile and megconnectomeRoot are required.');
end
if nargin < 3 || isempty(preStimTime), preStimTime = 1.5; end
if nargin < 4 || isempty(postStimTime), postStimTime = 2.5; end
if ~(ischar(rawFile) || isstring(rawFile)) || ~isfile(char(rawFile))
    error('MEG:EventRawMissing', 'Raw MEG file not found: %s', char(rawFile));
end
rawFile = char(rawFile);
megconnectomeRoot = char(megconnectomeRoot);
trialFunctionDir = fullfile(megconnectomeRoot, 'trial_functions');
if ~isfolder(trialFunctionDir)
    error('MEG:EventDecoderPath', 'Official trial_functions directory not found: %s', trialFunctionDir);
end
addpath(trialFunctionDir);

hdr = ft_read_header(rawFile);
recordingFsample = double(hdr.Fs);
recordingSamples = double(hdr.nSamples) * double(hdr.nTrials);
if ~isfinite(recordingFsample) || recordingFsample <= 0 || ...
        ~isfinite(recordingSamples) || recordingSamples < 1
    error('MEG:EventHeader', 'Invalid raw header for %s.', rawFile);
end

% Exact task block boundaries.  With zero padding, trl(:,1:2) are the
% official block start/end samples and columns 4:10 follow STIMBLOCKINFO.
cfg = make_cfg(rawFile, 'blocks', 0, 0);
[trlBlocks, ~, summaryBlocks, scanStart, scanEnd, warnBlocks] = trialfun_Wrkmem(cfg);

% Exact Fix15sec boundaries.  Columns 4:9 follow FIXATIONBLOCKINFO.
cfg = make_cfg(rawFile, 'fixation', 0, 0);
[trlFix, ~, summaryFix, scanStartFix, scanEndFix, warnFix] = trialfun_Wrkmem(cfg);

% Official image-locked task trials plus official Fix15sec pseudo-trials.
cfg = make_cfg(rawFile, 'trials', preStimTime, postStimTime);
[trlTrials, ~, summaryTrials, scanStartTrials, scanEndTrials, warnTrials] = trialfun_Wrkmem(cfg);

if any([scanStart scanStartFix scanStartTrials] ~= scanStart) || ...
        any([scanEnd scanEndFix scanEndTrials] ~= scanEnd)
    error('MEG:EventBoundsMismatch', 'Official trial modes returned inconsistent scan bounds.');
end
if scanStart < 1 || scanEnd > recordingSamples || scanEnd < scanStart
    error('MEG:EventBounds', 'Invalid official scan bounds [%g %g].', scanStart, scanEnd);
end

% Construct a common block record and restore temporal order.  Task blocks
% retain their official within-run block number; fixation blocks have no such
% number in the official fixation table and receive a fixation index.
blocks = struct('block_index', {}, 'block_number', {}, 'condition_code', {}, ...
    'image_type', {}, 'start_sample', {}, 'end_sample', {}, ...
    'previous_memory_type', {}, 'next_memory_type', {});
for i = 1:size(trlBlocks, 1)
    b.block_index = NaN;
    b.block_number = trlBlocks(i, 4);
    b.condition_code = trlBlocks(i, 6); % 1=0-back, 2=2-back
    b.image_type = trlBlocks(i, 5);
    b.start_sample = trlBlocks(i, 9);
    b.end_sample = trlBlocks(i, 10);
    b.previous_memory_type = trlBlocks(i, 8);
    b.next_memory_type = NaN;
    blocks(end+1) = b; %#ok<AGROW>
end
for i = 1:size(trlFix, 1)
    b.block_index = NaN;
    b.block_number = NaN;
    b.condition_code = 3; % Fix15sec
    b.image_type = 0;
    b.start_sample = trlFix(i, 8);
    b.end_sample = trlFix(i, 9);
    b.previous_memory_type = trlFix(i, 5);
    b.next_memory_type = trlFix(i, 7);
    blocks(end+1) = b; %#ok<AGROW>
end
if isempty(blocks)
    error('MEG:EventNoBlocks', 'Official decoder returned no task or fixation blocks.');
end
[~, order] = sort([blocks.start_sample]);
blocks = blocks(order);
for i = 1:numel(blocks)
    blocks(i).block_index = i;
end

% Look up parent block index from official block number or exact fixation
% boundary.  A trial is retained only when it has a valid parent block.
trials = struct('trial_index', {}, 'condition_code', {}, 'memory_type', {}, ...
    'image_type', {}, 'block_index', {}, 'block_number', {}, ...
    'trial_index_in_block', {}, 'event_sample', {}, 'window_start_sample', {}, ...
    'window_end_sample', {}, 'parent_start_sample', {}, 'parent_end_sample', {}, ...
    'is_pseudo_fixation', {});
preSamples = floor(recordingFsample * preStimTime);
for i = 1:size(trlTrials, 1)
    imageType = trlTrials(i, 7);
    memoryType = trlTrials(i, 8);
    isTask = isfinite(imageType) && imageType > 0 && ...
        isfinite(memoryType) && any(memoryType == [1 2]);
    isFix = isfinite(imageType) && imageType == 0;
    if ~(isTask || isFix)
        continue
    end
    if isTask
        blockNumber = trlTrials(i, 5);
        blockIndex = find([blocks.condition_code] == memoryType & ...
            [blocks.block_number] == blockNumber, 1, 'first');
        eventSample = trlTrials(i, 10);
    else
        blockNumber = NaN;
        fixStart = trlTrials(i, 10);
        fixEnd = trlTrials(i, 11);
        blockIndex = find([blocks.condition_code] == 3 & ...
            abs([blocks.start_sample] - fixStart) < 0.5 & ...
            abs([blocks.end_sample] - fixEnd) < 0.5, 1, 'first');
        eventSample = trlTrials(i, 1) + preSamples;
        memoryType = 3;
        imageType = 0;
    end
    if isempty(blockIndex)
        continue
    end
    tr.window_start_sample = trlTrials(i, 1);
    tr.window_end_sample = trlTrials(i, 2);
    tr.event_sample = eventSample;
    tr.condition_code = memoryType;
    tr.memory_type = memoryType;
    tr.image_type = imageType;
    tr.block_index = blockIndex;
    tr.block_number = blockNumber;
    tr.parent_start_sample = blocks(blockIndex).start_sample;
    tr.parent_end_sample = blocks(blockIndex).end_sample;
    tr.is_pseudo_fixation = isFix;
    tr.trial_index = numel(trials) + 1;
    tr.trial_index_in_block = NaN;
    trials(end+1) = tr; %#ok<AGROW>
end

% Assign a stable within-block index in event order.  The official task
% sequence is retained when available, but sorting is robust to the official
% task-then-fix output order.
for iBlock = 1:numel(blocks)
    ix = find([trials.block_index] == iBlock);
    if isempty(ix), continue; end
    [~, ord] = sort([trials(ix).event_sample]);
    ix = ix(ord);
    for j = 1:numel(ix)
        trials(ix(j)).trial_index_in_block = j;
    end
end

decoded = struct;
decoded.schema_version = 1;
decoded.raw_file = rawFile;
decoded.recording_fsample_hz = recordingFsample;
decoded.recording_samples = recordingSamples;
decoded.scan_start_sample = double(scanStart);
decoded.scan_end_sample = double(scanEnd);
decoded.blocks = blocks;
decoded.trials = trials;
decoded.trial_pre_seconds = preStimTime;
decoded.trial_post_seconds = postStimTime;
decoded.fixation_isi_seconds = 2.5;
decoded.summary_blocks = summaryBlocks;
decoded.summary_fixation = summaryFix;
decoded.summary_trials = summaryTrials;
decoded.warning_count = count_warnings({warnBlocks, warnFix, warnTrials});
end

function cfg = make_cfg(rawFile, cutMode, preStimTime, postStimTime)
cfg = [];
cfg.datafile = rawFile;
cfg.trialdef.lockedOn = 'Image';
cfg.trialdef.cutMode = cutMode;
cfg.trialdef.preStimTime = preStimTime;
cfg.trialdef.postStimTime = postStimTime;
end

function n = count_warnings(warnings)
n = 0;
for i = 1:numel(warnings)
    w = warnings{i};
    if isempty(w), continue; end
    if iscell(w), n = n + numel(w); else, n = n + 1; end
end
end
