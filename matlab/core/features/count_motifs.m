function [numpattern1, numpattern2] = count_motifs(Map,size_threshold)
if nargin < 2
    size_threshold = 10;
end
ind = find(Map==1);
if length(ind)<size_threshold
    numpattern1=0;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern1 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern1);
    ps_k=zeros(numpattern1,1);
    for k = numpattern1:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    numpattern1 = numpattern1 - length(find(ps_k<size_threshold));
end

ind = find(Map==-1);
if length(ind)<size_threshold
    numpattern2=0;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern2 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern2);
    ps_k=zeros(numpattern2,1);
    for k = numpattern2:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    numpattern2 = numpattern2 - length(find(ps_k<size_threshold));
end

end