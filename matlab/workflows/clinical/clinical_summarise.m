function S=clinical_summarise(D)
%CLINICAL_SUMMARISE Published group comparisons, joint score and quadrants.
% Ordinary covariate-adjusted OLS follows the original figure producer.
% Descriptive means/SEM and Cohen's d use all finite metric observations;
% the adjusted model additionally requires complete covariates.
required={'PQ','Acc','group01','covMat','conditionNames','subjectList'};
assert(all(isfield(D,required)),'clinical:Fields','Incomplete cohort structure.');
n=numel(D.group01); cnum=numel(D.conditionNames);
assert(isequal(size(D.PQ),[n cnum]) && isequal(size(D.Acc),[n cnum]) && size(D.covMat,1)==n,'clinical:Shape','Cohort dimensions disagree.');
assert(numel(D.subjectList)==n && numel(unique(string(D.subjectList)))==n,'clinical:IDs','Participant IDs must be unique.');
g=D.group01(:);assert(all(ismember(g,[0 1])),'clinical:Group','Group must be HC=0, SCZ=1.');
[q,cutoffs]=clinical_quadrants(D.PQ,D.Acc,g);
effects=table(); proportions=table(); thresholds=table();
joint=nan(size(D.PQ));
for c=1:cnum
    joint(:,c)=finite_zscore(D.PQ(:,c))+finite_zscore(D.Acc(:,c));
    values=[D.PQ(:,c),D.Acc(:,c),joint(:,c)];
    names=["PQ","accuracy","joint_score"];
    for j=1:3
        y=values(:,j); hc=y(g==0 & isfinite(y)); scz=y(g==1 & isfinite(y));
        nh=numel(hc);ns=numel(scz); delta=mean(scz)-mean(hc);
        d=nan;
        if nh>1 && ns>1
            pooled=sqrt(((nh-1)*var(hc)+(ns-1)*var(scz))/(nh+ns-2));
            if pooled>eps,d=delta/pooled;end
        end
        valid=isfinite(y)&all(isfinite(D.covMat),2);
        x=[ones(sum(valid),1),g(valid),D.covMat(valid,:)];
        beta=nan;p=nan; status="insufficient_complete_cases";
        if sum(valid)>=8 && numel(unique(g(valid)))==2
            assert(rank(x)==size(x,2),'clinical:Rank','Rank-deficient adjusted group model.');
            assert(size(x,1)>size(x,2),'clinical:DF','No residual degrees of freedom.');
            mdl=fitlm(x(:,2:end),y(valid));
            beta=mdl.Coefficients.Estimate(2);p=mdl.Coefficients.pValue(2);status="ok";
        end
        row=table(string(D.conditionNames(c)),names(j),nh,ns,mean(hc),mean(scz),std(hc)/sqrt(nh),std(scz)/sqrt(ns),delta,d,sum(valid),beta,p,status,...
            'VariableNames',{'condition','metric','n_hc','n_scz','mean_hc','mean_scz','sem_hc','sem_scz','delta_scz_minus_hc','cohens_d','n_model','beta_group','p_group','status'});
        effects=[effects;row]; %#ok<AGROW>
    end
    thresholds=[thresholds;table(string(D.conditionNames(c)),cutoffs(c,1),cutoffs(c,2),'VariableNames',{'condition','pq_hc_p10','accuracy_hc_p10'})]; %#ok<AGROW>
    for group=0:1
        count=arrayfun(@(k)sum(q(:,c)==k & g==group),1:4);
        denominator=sum(count);
        for k=1:4
            labels=["neither_low","low_pq_only","low_accuracy_only","both_low"];
            proportions=[proportions;table(string(D.conditionNames(c)),group,labels(k),count(k),denominator,count(k)/denominator,...
                'VariableNames',{'condition','group01','quadrant','count','n_complete','proportion'})]; %#ok<AGROW>
        end
    end
end
% Preserve the original within-dataset family (3 metrics x all task loads).
p=effects.p_group; adjusted=nan(size(p)); ii=find(isfinite(p));
[ordered,order]=sort(p(ii)); m=numel(ordered);
if m>0,adjusted(ii(order))=min(1,flipud(cummin(flipud(ordered.*m./(1:m)'))));end
effects.q_group=adjusted;
S=struct('group_effects',effects,'thresholds',thresholds,'proportions',proportions,'quadrant',q,'joint_score',joint);
end

function z=finite_zscore(x)
z=nan(size(x));valid=isfinite(x);
if sum(valid)<2,return;end
sd=std(x(valid));if ~isfinite(sd)||sd<=eps,return;end
z(valid)=(x(valid)-mean(x(valid)))/sd;
end
