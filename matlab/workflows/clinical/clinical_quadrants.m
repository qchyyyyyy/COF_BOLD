function [assignment,thresholds]=clinical_quadrants(pq,accuracy,group)
%CLINICAL_QUADRANTS HC tenth-percentile partition used in Fig. 5 and ED11.
% Columns are conditions; rows are participants. Missing pairs remain NaN.
% 1=neither low, 2=low PQ only, 3=low accuracy only, 4=both low.
assert(isequal(size(pq),size(accuracy)) && size(pq,1)==numel(group),'clinical:Shape','Input dimensions disagree.');
assert(all(ismember(group,[0 1])),'clinical:Group','Group must be HC=0, SCZ=1.');
group=group(:); assignment=nan(size(pq)); thresholds=nan(size(pq,2),2);
for c=1:size(pq,2)
    valid=isfinite(pq(:,c)) & isfinite(accuracy(:,c));
    hc=valid & group==0;
    assert(any(hc),'clinical:HC','No complete HC pairs for condition %d.',c);
    thresholds(c,:)=[prctile(pq(hc,c),10),prctile(accuracy(hc,c),10)];
    assignment(valid,c)=1+double(pq(valid,c)<thresholds(c,1))+2*double(accuracy(valid,c)<thresholds(c,2));
end
end
