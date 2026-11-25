function [sizepattern1, sizepattern2] = measure_motif_sizes(Map)
ind = find(Map==1);

if isempty(ind)
    sizepattern1 = 0;
elseif length(ind)==1
    sizepattern1 = 1;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern1 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern1);
    ps_k=zeros(numpattern1,1);
    for k = numpattern1:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    sizepattern1 = mean(ps_k);
end

ind = find(Map==-1);

if isempty(ind)
    sizepattern2 = 0;
elseif length(ind)==1
    sizepattern2 = 1;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern2 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern2);
    ps_k=zeros(numpattern2,1);
    for k = numpattern2:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    sizepattern2 = mean(ps_k);
end

end