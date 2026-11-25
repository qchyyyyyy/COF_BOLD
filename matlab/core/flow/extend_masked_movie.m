function [extendedMovie, extendedMask] = extend_masked_movie(movie, mask, nRings)
%HARMONIC_EXTEND_MASKED_MOVIE Fill exterior mask rings by a Laplace solve.

validateattributes(movie, {'single', 'double'}, {'nonempty'}, ...
    mfilename, 'movie');
validateattributes(mask, {'numeric', 'logical'}, ...
    {'2d', 'nonempty', 'finite'}, mfilename, 'mask');
validateattributes(nRings, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'nRings');
if size(movie, 1) ~= size(mask, 1) || size(movie, 2) ~= size(mask, 2)
    error('COF:HarmonicExtensionSizeMismatch', ...
        'The first two movie dimensions must match the mask size.');
end

mask = logical(mask);
if ~any(mask(:))
    error('COF:HarmonicExtensionEmptyMask', ...
        'The original mask must not be empty.');
end

[nRows, nCols, nFrames] = size(movie);
movieMatrix = reshape(movie, nRows * nCols, nFrames);
insideValues = movieMatrix(mask(:), :);
if any(~isfinite(insideValues(:)))
    error('COF:HarmonicExtensionNonfiniteInput', ...
        'Movie values inside the original mask must be finite.');
end

extendedMask = mask;
for iRing = 1:nRings
    extendedMask = conv2(double(extendedMask), ones(3), 'same') > 0;
end
extendedMovie = movie;
ringMask = extendedMask & ~mask;
ringIndices = find(ringMask);
if isempty(ringIndices)
    return
end

nUnknown = numel(ringIndices);
unknownLabels = zeros(nRows, nCols);
unknownLabels(ringIndices) = 1:nUnknown;
rowIndices = zeros(5 * nUnknown, 1);
columnIndices = zeros(5 * nUnknown, 1);
matrixValues = zeros(5 * nUnknown, 1);
rhs = zeros(nUnknown, nFrames);
if ~isreal(movie)
    rhs = complex(rhs);
end
nEntries = 0;
offsets = [-1 0; 1 0; 0 -1; 0 1];

for iUnknown = 1:nUnknown
    [row, col] = ind2sub([nRows, nCols], ringIndices(iUnknown));
    degree = 0;
    for iOffset = 1:size(offsets, 1)
        neighborRow = row + offsets(iOffset, 1);
        neighborCol = col + offsets(iOffset, 2);
        inImage = neighborRow >= 1 && neighborRow <= nRows && ...
            neighborCol >= 1 && neighborCol <= nCols;
        if ~inImage || ~extendedMask(neighborRow, neighborCol)
            continue
        end

        degree = degree + 1;
        neighborIndex = sub2ind([nRows, nCols], neighborRow, neighborCol);
        neighborLabel = unknownLabels(neighborIndex);
        if neighborLabel > 0
            nEntries = nEntries + 1;
            rowIndices(nEntries) = iUnknown;
            columnIndices(nEntries) = neighborLabel;
            matrixValues(nEntries) = -1;
        else
            rhs(iUnknown, :) = rhs(iUnknown, :) + ...
                double(movieMatrix(neighborIndex, :));
        end
    end

    nEntries = nEntries + 1;
    rowIndices(nEntries) = iUnknown;
    columnIndices(nEntries) = iUnknown;
    matrixValues(nEntries) = degree;
end

laplacian = sparse(rowIndices(1:nEntries), columnIndices(1:nEntries), ...
    matrixValues(1:nEntries), nUnknown, nUnknown);
ringValues = laplacian \ rhs;
movieMatrix(ringIndices, :) = cast(ringValues, 'like', movie);
extendedMovie = reshape(movieMatrix, size(movie));
end
