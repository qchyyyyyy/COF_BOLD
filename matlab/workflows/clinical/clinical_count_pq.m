function counts=clinical_count_pq(left,right,mask,method)
%CLINICAL_COUNT_PQ Original global source/sink counts from bilateral VVFs.
% Arrays are grid-y x grid-x x time. For PG, preserve historical half-grid
% concatenation followed by 2 x 2 x 2 valid averaging before differentiation.
method=validatestring(method,{'COF','OF','PG'});
% Multiple runs are joined BEFORE pooled derivative standardization, as in
% the original ds000115 producer. PG averaging is performed within each run.
if ~iscell(left),left={left};end
if ~iscell(right),right={right};end
assert(numel(left)==numel(right)&&~isempty(left),'clinical:Runs','Hemisphere run lists disagree.');
ux=[];uy=[];
for k=1:numel(left)
    l=left{k};r=right{k};
    assert(all(isfield(l,{'Ux','Uy'}))&&all(isfield(r,{'Ux','Uy'})),'clinical:Velocity','Missing velocity components.');
    assert(isequal(size(l.Ux),size(l.Uy))&&isequal(size(r.Ux),size(r.Uy)),'clinical:Shape','Velocity components disagree.');
    if strcmp(method,'PG')
        x=cat(2,l.Ux(:,1:end-1,:),r.Ux);y=cat(2,l.Uy(:,1:end-1,:),r.Uy);
        x=convn(x,ones(2,2,2)/8,'valid');y=convn(y,ones(2,2,2)/8,'valid');
    else
        x=cat(2,l.Ux,r.Ux);y=cat(2,l.Uy,r.Uy);
    end
    ux=cat(3,ux,x);uy=cat(3,uy,y);
end
assert(isequal(size(ux,[1 2]),size(mask)),'clinical:Shape','Bilateral velocity grid does not match mask.');
[div,~]=velocity_derivatives(ux,uy,mask,3,1,1,2);
count=size(div,3); source=zeros(1,count);sink=zeros(1,count);
for t=1:count,[source(t),sink(t)]=count_motifs(div(:,:,t),10);end
counts=struct('NumSource',source,'NumSink',sink,'PQ',source+sink);
end
