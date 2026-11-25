function S=clinical_disorganization_trend(y,quadrant,group,covariates)
%CLINICAL_DISORGANIZATION_TREND Original Fig. 5g two-sided ordinary OLS trend.
% y is z_dis_4grp; only SCZ participants are fitted. Covariates are age,
% education and binary sex. quadrant is the FROZEN historical Fig. 5g order
% (1 neither low, 2 low PQ only, 3 low accuracy only, 4 both low).
% No other phenotype tests, omnibus model, robust SE or BH family is run.
y=y(:);quadrant=quadrant(:);group=group(:);
assert(numel(y)==numel(quadrant)&&numel(y)==numel(group)&&isequal(size(covariates),[numel(y),3]),'clinical:Shape','Trend input dimensions disagree.');
assert(all(isnan(quadrant)|ismember(quadrant,1:4)) && all(ismember(group,[0 1])),'clinical:Labels','Invalid quadrant/group labels.');
valid=group==1 & isfinite(y) & isfinite(quadrant) & all(isfinite(covariates),2);
yy=y(valid); rr=quadrant(valid);
assert(numel(yy)>=6 && numel(unique(rr))>=2,'clinical:Sample','Insufficient data for trend.');
yy=standardize(yy);rr=standardize(rr);
x=[ones(sum(valid),1),rr,covariates(valid,:)];
assert(rank(x)==size(x,2),'clinical:Rank','Rank-deficient historical trend model.');
assert(size(x,1)>size(x,2),'clinical:DF','No residual degrees of freedom.');
mdl=fitlm(x(:,2:end),yy);
S=table(sum(valid),mdl.Coefficients.Estimate(2),mdl.Coefficients.pValue(2),mdl.DFE,...
    'VariableNames',{'n','beta_trend','p_trend','residual_df'});
for k=1:4
    v=y(valid & quadrant==k);
    S.(sprintf('n_quadrant%d',k))=numel(v);
    S.(sprintf('mean_quadrant%d',k))=mean(v);
    S.(sprintf('sem_quadrant%d',k))=std(v)/sqrt(numel(v));
end
end

function z=standardize(x)
sd=std(x);if sd<=eps || ~isfinite(sd),sd=1;end
z=(x-mean(x))/sd;
end
