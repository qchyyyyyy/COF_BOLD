function test_core_clinical()
%TEST_CORE_CLINICAL Synthetic numerical contracts; no participant data needed.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup('clinical');
rng(20261125,'twister');
assert(abs(wrap_angle_difference(pi,0)+pi)<1e-12);
[xx,yy]=meshgrid(1:17);phase=.1*xx+.05*yy;
[gx,gy]=phase_gradient(repmat(phase,1,1,3));
assert(max(abs(gx(2:16,2:16,:)+.1),[],'all')<1e-12);
assert(max(abs(gy(2:16,2:16,:)+.05),[],'all')<1e-12);
p=struct('rho',5e-4,'tau',1e-8,'sigma',8);
[vx,vy]=cof_solve({repmat(cos(phase),1,1,3),repmat(sin(phase),1,1,3)},true(17),p);
assert(isequal(size(vx),[16 16 2])&&isequal(size(vy),size(vx)));
assert(all(isfinite(vx),'all')&&all(isfinite(vy),'all'));
assert(max(abs(vx),[],'all')<1e-10&&max(abs(vy),[],'all')<1e-10);
% Previously broken optional unsmoothed derivative branch now has a contract.
[dd,cc]=velocity_derivatives(xx,yy,true(17));
assert(max(abs(dd-2),[],'all')<1e-12);assert(isequal(size(cc),[17 17]));
map=zeros(10);map(1:3,1:4)=1;map(6:8,6:9)=-1;
[ns,nk]=count_motifs(map);assert(ns==1&&nk==1);
% Strict boundary and missing-pair handling.
pq=[(1:20)';1.5;1;5;1;NaN];acc=[(1:20)';1.5;5;1;1;4];group=[zeros(20,1);ones(5,1)];
[q,cut]=clinical_quadrants(pq,acc,group);
assert(isequal(cut,[2.5 2.5]));assert(isequal(q(22:24),[2;3;4])&&isnan(q(25)));
[qq,~]=clinical_quadrants([pq;2.5],[acc;2.5],[group;1]);assert(qq(end)==1);
n=80;g=[zeros(n/2,1);ones(n/2,1)];cov=[20+20*rand(n,1),rand(n,1),double(rand(n,1)>.5)];
D=struct('PQ',randn(n,2)+g*.4,'Acc',randn(n,2)-g*.5,'group01',g,'covMat',cov,'conditionNames',["0-back","2-back"],'subjectList',"synthetic_"+string((1:n)'));
S=clinical_summarise(D);assert(height(S.group_effects)==6&&height(S.proportions)==16);
assert(max(abs(S.joint_score(:,1)-(zscore(D.PQ(:,1))+zscore(D.Acc(:,1)))))<1e-12);
model=fitlm([g cov],D.PQ(:,1));
assert(abs(S.group_effects.beta_group(1)-model.Coefficients.Estimate(2))<1e-12);
assert(abs(S.group_effects.p_group(1)-model.Coefficients.pValue(2))<1e-12);
assert(all(S.group_effects.q_group>=S.group_effects.p_group-1e-12));
% Independent closed-form OLS comparison for Fig. 5g.
quadrant=repmat((1:4)',n/4,1);y=.2*quadrant+.1*cov(:,1)+randn(n,1);
trend=clinical_disorganization_trend(y,quadrant,g,cov);
v=g==1;x=[ones(sum(v),1),zscore(quadrant(v)),cov(v,:)];z=zscore(y(v));
b=x\z;df=size(x,1)-size(x,2);e=z-x*b;se=sqrt((e'*e/df)*inv(x'*x));
pr=2*tcdf(-abs(b(2)/se(2,2)),df);
assert(abs(trend.beta_trend-b(2))<1e-10&&abs(trend.p_trend-pr)<1e-10);
assert(~ismember('p_regime',trend.Properties.VariableNames)&&~ismember('q_trend',trend.Properties.VariableNames));
% Runtime input guards must fail before writing any partial subject data.
expect_error(@()clinical_quadrants(pq,acc,ones(size(group))),'clinical:HC');
expect_error(@()clinical_disorganization_trend(y,quadrant,g,zeros(n,3)),'clinical:Rank');
expect_error(@()clinical_output_directory(root),'clinical:Output');
fprintf('PASS core/clinical: phase, stationary COF, derivatives, motifs, strict/missing quadrants, joint score, OLS, BH, trend and input guards.\n');
end

function expect_error(fn,id)
try,fn();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.identifier);return;end
error('test:MissingError','Expected %s.',id);
end
