function S=clinical_load_comparison(pq,events,group)
%CLINICAL_LOAD_COMPARISON Descriptive global PQ by condition (ED12).
% pq and events: participants x time, aligned after the historical exclusion.
% Call once per OF, PG and COF input; events 0/1/2=fixation/0-back/2-back.
assert(isequal(size(pq),size(events))&&size(pq,1)==numel(group),'clinical:Shape','Input dimensions disagree.');
assert(all(ismember(group,[0 1])),'clinical:Group','Group must be HC=0, SCZ=1.');
means=nan(size(pq,1),3);S=table();names=["Fixation","0-back","2-back"];
for c=0:2
    x=pq;x(events~=c)=nan;means(:,c+1)=mean(x,2,'omitnan');
    for g=0:1
        v=means(group(:)==g,c+1);v=v(isfinite(v));
        S=[S;table(names(c+1),g,numel(v),mean(v),std(v)/sqrt(numel(v)),...
            'VariableNames',{'condition','group01','n','mean_pq','sem_pq'})]; %#ok<AGROW>
    end
end
end
