function [Patterns, DwellTime] = track_motifs(Map,value,size_threshold)
if nargin < 3
    size_threshold = 10;
end
Cluster = zeros(size(Map));
T=size(Cluster,3);
for t = 1:T
    ind = find(Map(:,:,t)==value);
    if length(ind)>=size_threshold
        [I,J] = ind2sub(size(Map(:,:,t)),ind);
        Z = linkage([I J],'single');
        K = length(find(Z(:,3)>sqrt(2)))+1;
        Cluster_ind = cluster(Z,'MaxClust',K);
        Cluster_map = zeros(size(Cluster(:,:,t)));
        Cluster_map(ind) = Cluster_ind;
        Cluster(:,:,t) = Cluster_map;
    end
end
clearvars Map Cluster_map Cluster_ind Z K I J


visited = cell([T 1]);
for t = 1:T
    K = max(Cluster(:,:,t),[],'all');
    visited{t} = zeros([K 1]);
end
Patterns.ind=[];
DwellTime = [];
cnt = 1;
for t = 1:T-1
    K = max(Cluster(:,:,t),[],'all');
    for cc = 1:K
        if ~visited{t}(cc)
            ind_t = find(Cluster(:,:,t)==cc);
            Ind{1} = ind_t;
            [Patterns.ind{cnt},visited] = match_cluster(Ind,cc,t,visited,Cluster);
            DwellTime(cnt) = length(Patterns.ind{cnt});
            Patterns.time{cnt} = [t t+length(Patterns.ind{cnt})-1];
            cnt = cnt + 1;
        end
    end
end
clearvars Cluster
end % end function

function [Ind,visited] = match_cluster(Ind,cluster,t,visited,Cluster)
T = size(Cluster,3);
visited{t}(cluster)=1;
if t == T
    return;
end
for nextcluster = 1:max(Cluster(:,:,t+1),[],'all')
    if ~visited{t+1}(nextcluster)
        ind_nextcluster = find(Cluster(:,:,t+1)==nextcluster);
        if ~isempty(intersect(Ind{end},ind_nextcluster))
            Ind{length(Ind)+1} = ind_nextcluster;
            [Ind,visited] = match_cluster(Ind,nextcluster,t+1,visited,Cluster);
            break;
        end
    end
end

end