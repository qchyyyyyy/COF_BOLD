function [sizepattern1, sizepattern2] = meg_measure_motif_sizes(Map)
%MEG_MEASURE_MOTIF_SIZES Mean positive and negative pattern size in pixels.
sizepattern1 = mean_sign(Map, 1);
sizepattern2 = mean_sign(Map, -1);
end

function value = mean_sign(Map, signValue)
ind = find(Map == signValue);
if isempty(ind)
    value = 0;
elseif numel(ind) == 1
    value = 1;
else
    [I,J] = ind2sub(size(Map), ind);
    Z = linkage([I J], 'single');
    nPatterns = nnz(Z(:,3) > sqrt(2)) + 1;
    clusterIndex = cluster(Z, 'MaxClust', nPatterns);
    clusterSize = zeros(nPatterns,1);
    for k = nPatterns:-1:1
        clusterSize(k) = nnz(clusterIndex == k);
    end
    value = mean(clusterSize);
end
end
