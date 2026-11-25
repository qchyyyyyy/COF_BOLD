function test_clinical_inputs()
%TEST_CLINICAL_INPUTS Across-run normalisation and participant alignment guards.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));addpath(root);setup('clinical');
tmp=tempname;mkdir(tmp);cleanup=onCleanup(@()rmdir(tmp,'s')); %#ok<NASGU>
for method=["COF","PG"]
    if method=="PG",ny=29;nx=22;nt=4;else,ny=28;nx=21;nt=3;end
    [xx,yy]=meshgrid((1:nx)-mean(1:nx),(1:ny)-mean(1:ny));
    envelope=exp(-(xx.^2+yy.^2)/50);
    vx=xx.*envelope;vy=yy.*envelope;
    scales=[.01 1 20];left=cell(1,3);right=left;ux=[];uy=[];
    separate=[];
    mask=ones(28,42);
    for k=1:3
        field=struct('Ux',scales(k)*repmat(vx,1,1,nt),'Uy',scales(k)*repmat(vy,1,1,nt));
        left{k}=field;right{k}=field;
        if method=="PG"
            x=convn(cat(2,field.Ux(:,1:end-1,:),field.Ux),ones(2,2,2)/8,'valid');
            y=convn(cat(2,field.Uy(:,1:end-1,:),field.Uy),ones(2,2,2)/8,'valid');
        else
            x=cat(2,field.Ux,field.Ux);y=cat(2,field.Uy,field.Uy);
        end
        ux=cat(3,ux,x);uy=cat(3,uy,y);
        singleRun=clinical_count_pq(field,field,mask,method);
        separate=[separate singleRun.PQ]; %#ok<AGROW>
    end
    actual=clinical_count_pq(left,right,mask,method);
    % Historical producer concatenates prepared velocity runs before calc_DivCurl.
    [d,~]=velocity_derivatives(ux,uy,mask,3,1,1,2);
    ns=zeros(1,size(d,3));nk=ns;
    for t=1:size(d,3),[ns(t),nk(t)]=count_motifs(d(:,:,t),10);end
    assert(isequal(actual.NumSource,ns)&&isequal(actual.NumSink,nk));
    assert(numel(actual.PQ)==9,'test:RunBoundary','PG averaging crossed a run boundary.');
    assert(~isequal(actual.PQ,separate),'test:WeakFixture','Fixture must distinguish pooled and per-run thresholds.');
    % Exercise the documented nested JSON run matrix through file aggregation.
    cii_mask_left=ones(29,22);cii_mask_right=cii_mask_left; %#ok<NASGU>
    save(fullfile(tmp,'coordinate_cii_v2.mat'),'cii_mask_left','cii_mask_right');
    dirs=cell(1,3);
    for k=1:3
        dirs{k}=fullfile(tmp,sprintf('%s_run%d',method,k));mkdir(dirs{k});
        VelField=left{k};save(fullfile(dirs{k},sprintf('VelField_lh_%s.mat',method)),'VelField');
        VelField=right{k};save(fullfile(dirs{k},sprintf('VelField_rh_%s.mat',method)),'VelField');
    end
    row=['[' strjoin(cellfun(@jsonencode,dirs,'UniformOutput',false),',') ']'];
    out=fullfile(tmp,['aggregate_' char(method)]);
    payload=sprintf('{"subject_ids":["testA","testB"],"run_dirs":[%s,%s],"frames_per_run":3,"method":%s,"template_root":%s,"output_dir":%s,"output_name":"counts.mat"}', ...
        row,row,jsonencode(char(method)),jsonencode(tmp),jsonencode(out));
    configFile=fullfile(tmp,['aggregate_' char(method) '.json']);
    fid=fopen(configFile,'w');fwrite(fid,payload,'char');fclose(fid);
    clinical_aggregate_pq(configFile);
    combined=load(fullfile(out,'counts.mat'));
    assert(isequal(string(combined.subjectList),["testA";"testB"]));
    assert(isequal(combined.NumSource,repmat(ns,2,1))&&isequal(combined.NumSink,repmat(nk,2,1)));
end

% Minimal ds000115-style inputs, including one excluded row.
resultsRoot=fullfile(tmp,'results');site=fullfile(resultsRoot,'SCZ_WM');mkdir(site);mkdir(fullfile(site,'group'));
participant_id=["test01";"test02";"test03";"test04"];
condit=["CON";"SCZ";"CON";"SCZ"];age=[22;31;27;35];yrschool=[14;12;16;15];gender=["MALE";"FEMALE";"FEMALE";"MALE"];
T=table(participant_id,condit,age,yrschool,gender);
for field=["nback0_targ","nback0_nont","nback1_targ","nback1_nont","nback2_targ","nback2_nont"]
    T.(field)=[.9;.7;.8;.6];
end
writetable(T,fullfile(site,'participants.tsv'),'FileType','text','Delimiter','\t');
cfg=struct('dataRoot',tmp,'resultsRoot',resultsRoot,'method','COF','hemisphere','LR', ...
    'excludedRows',struct('ZIB',[],'ds000115',2),'allow_historical_positional_input',false);
patternFile=fullfile(site,'group','numPattern_global_LR_COF.mat');
NumSource=repmat((1:4)',1,1425);NumSink=zeros(4,1425);subjectList=participant_id;
save(patternFile,'NumSource','NumSink','subjectList');
D=clinical_collect_cohort(cfg,'ds000115');
assert(isequal(D.subjectList,participant_id([1 3 4]))&&isequal(D.PQ,repmat([1;3;4],1,3)));
subjectList=subjectList([3 2 1 4]);save(patternFile,'NumSource','NumSink','subjectList');
expect_error(@()clinical_collect_cohort(cfg,'ds000115'),'clinical:Alignment');
subjectList=participant_id([1 3 4]);NumSource=NumSource([1 3 4],:);NumSink=NumSink([1 3 4],:);
save(patternFile,'NumSource','NumSink','subjectList');
D=clinical_collect_cohort(cfg,'ds000115');assert(isequal(D.PQ,repmat([1;3;4],1,3)));
save(patternFile,'NumSource','NumSink');
expect_error(@()clinical_collect_cohort(cfg,'ds000115'),'clinical:Alignment');
cfg.allow_historical_positional_input=true;
D=clinical_collect_cohort(cfg,'ds000115');assert(isequal(D.PQ,repmat([1;3;4],1,3)));
fprintf('PASS clinical inputs: pooled thresholds, within-run PG averaging, historical exclusion, reordered IDs and explicit positional fallback.\n');
end

function expect_error(fn,id)
try,fn();catch ME,assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.identifier);return;end
error('test:MissingError','Expected %s.',id);
end
