function [power, validMask, metadata, manifest] = meg_load_power_session(sourceDir)
%MEG_LOAD_POWER_SESSION Assemble nonoverlapping chunk cores in memory.

manifestFile = fullfile(sourceDir, 'session_manifest.json');
metadataFile = fullfile(sourceDir, 'session_metadata.mat');
if ~isfile(manifestFile) || ~isfile(metadataFile)
    error('MEG:PowerSessionFiles', 'Missing manifest or metadata in %s.', sourceDir);
end
manifest = jsondecode(fileread(manifestFile));
loaded = load(metadataFile, 'metadata');
metadata = loaded.metadata;
nSources = double(manifest.source_count);
nFrames = double(manifest.total_frames);
power = nan(nSources, nFrames, 'single');
validMask = false(1, nFrames);

chunks = manifest.chunks;
for iChunk = 1:numel(chunks)
    chunk = load(fullfile(sourceDir, chunks(iChunk).filename), ...
        'source_power', 'sample_index', 'core_time_mask', 'valid_time_mask');
    corePosition = find(chunk.core_time_mask);
    frameIndex = double(chunk.sample_index(corePosition));
    if any(any(isfinite(power(:, frameIndex))))
        error('MEG:ChunkCoreOverlap', 'Chunk cores overlap in %s.', sourceDir);
    end
    power(:, frameIndex) = chunk.source_power(:, corePosition);
    validMask(frameIndex) = chunk.valid_time_mask(corePosition);
end
if nnz(validMask) ~= double(manifest.valid_frames)
    error('MEG:ValidFrameCount', 'Manifest and assembled valid-frame counts disagree.');
end
end
