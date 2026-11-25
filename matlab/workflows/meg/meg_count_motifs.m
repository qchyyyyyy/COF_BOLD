function [numpattern1, numpattern2] = meg_count_motifs(Map, size_threshold)
%MEG_COUNT_MOTIFS Count positive and negative connected spatial patterns.
% Project 01 convention: single-linkage distance threshold sqrt(2), with a
% default minimum pattern size of 10 pixels.
if nargin < 2, size_threshold = 10; end
numpattern1 = count_sign(Map, 1, size_threshold);
numpattern2 = count_sign(Map, -1, size_threshold);
end

function count = count_sign(Map, signValue, size_threshold)
ind = find(Map == signValue);
if numel(ind) < size_threshold
    count = 0;
    return
end
[I,J] = ind2sub(size(Map), ind);
Z = linkage([I J], 'single');
count = nnz(Z(:,3) > sqrt(2)) + 1;
clusterIndex = cluster(Z, 'MaxClust', count);
clusterSize = zeros(count,1);
for k = count:-1:1
    clusterSize(k) = nnz(clusterIndex == k);
end
count = count - nnz(clusterSize < size_threshold);
end
