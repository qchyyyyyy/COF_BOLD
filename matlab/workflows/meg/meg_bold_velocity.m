function summary = meg_bold_velocity(configFile,subjectId,scanId,mode)
%MEG_BOLD_VELOCITY Hemisphere-first gap-aware COF for BOLD-like signals.

if nargin < 4, mode='locked'; end
cfg = meg_bold_config(configFile);
subjectId=char(subjectId); scanId=char(scanId); mode=char(mode);
assert(strcmp(mode,'locked'),'MEG:Parameters','Only the manuscript parameter set is supported.');
assert(exist('nanconv','file')==2,'Run setup(''meg'') first.');
preprocFile=fullfile(cfg.processed_root,subjectId,'Wrkmem',char(cfg.derived_subdir), ...
    ['scan-' scanId],'preprocessed_signal.mat');
if ~isfile(preprocFile), error('MEG:WorkflowPreprocMissing','Missing %s.',preprocFile); end
loaded=load(preprocFile,'data'); data=loaded.data;
coordinates=load(cfg.coordinate_file,'cii_mask_left','cii_mask_right');
parameterSets=build_parameters(cfg,mode);
branches=normalise_branches(cfg.cof.branches);
rows=struct('subject_id',{},'scan_id',{},'branch',{},'parameter_name',{}, ...
    'rho',{},'tau',{},'sigma',{},'finite_fraction',{},'speed_median',{}, ...
    'speed_p99',{},'spatial_roughness',{},'brightness_residual_ratio',{}, ...
    'intervals',{},'output_directory',{});
for iBranch=1:numel(branches)
    branch=branches{iBranch};
    brainL=data.(['BrainImg_lh_' branch]);
    brainR=data.(['BrainImg_rh_' branch]);
    for iParam=1:numel(parameterSets)
        p=parameterSets(iParam);
        outputDir=fullfile(fileparts(preprocFile),['cof_' branch '_' p.name]);
        if isfile(fullfile(outputDir,'cof_manifest.json')) && cfg.allow_existing_output
            manifest=jsondecode(fileread(fullfile(outputDir,'cof_manifest.json')));
            validate_existing_manifest(manifest,branch,p,data);
            diag=manifest.diagnostics;
        else
            if isfolder(outputDir)
                contents=dir(outputDir); contents=contents(~ismember({contents.name},{'.','..'}));
                if ~isempty(contents), error('MEG:WorkflowExistingCOF','Refusing partial output %s.',outputDir); end
            else, mkdir(outputDir); end
            [manifest,diag]=compute_cof(brainL,brainR,coordinates,data,p,cfg,outputDir,branch);
        end
        rows(end+1)=make_row(subjectId,scanId,branch,p,diag,outputDir); %#ok<AGROW>
        if strcmp(mode,'locked')
            metricDir=fullfile(outputDir,'metrics');
            if ~isfile(fullfile(metricDir,'divcurl_manifest.json'))
                if isfolder(metricDir)
                    error('MEG:WorkflowPartialMetrics', ...
                        'Metric directory exists without a complete manifest: %s',metricDir);
                end
                options=struct('filtSigma',cfg.metrics.filt_sigma, ...
                    'doZscore',logical(cfg.metrics.do_zscore), ...
                    'doBinarize',logical(cfg.metrics.do_binarize), ...
                    'nStd',cfg.metrics.n_std,'sizeThreshold',cfg.metrics.size_threshold, ...
                    'saveMaps',logical(cfg.metrics.save_maps_for_cross_modal));
                meg_native_motif_metrics(outputDir,metricDir,[1 Inf],options);
            end
            validate_metric_manifest(metricDir,p,cfg);
            meg_bold_export_metrics(configFile,subjectId,scanId,branch,p.name,outputDir,metricDir);
        end
    end
end
if strcmp(mode,'locked')
    if numel(parameterSets)~=1
        error('MEG:WorkflowLockedCardinality','Locked mode must contain exactly one parameter set.');
    end
    parameterId=parameterSets(1).name;
    resultDir=cof_result_dir(cfg,parameterId,subjectId,scanId);
end
if ~isfolder(resultDir),mkdir(resultDir);end
diagnostics=struct2table(rows);
writetable(diagnostics,fullfile(resultDir,['cof_' mode '_diagnostics.csv']));
summary=struct('schema_version',1,'status','complete','workflow','meg_bold', ...
    'subject_id',subjectId,'scan_id',scanId,'mode',mode,'parameter_id',parameterId, ...
    'branches',{branches}, ...
    'parameter_count',numel(parameterSets),'diagnostic_rows',height(diagnostics), ...
    'created',datestr(now,30));
write_json(fullfile(resultDir,['cof_' mode '_manifest.json']),summary);
end

function parameters=build_parameters(cfg,mode)
assert(strcmp(mode,'locked'),'MEG:Parameters','Only locked parameters are supported.');
q=cfg.cof.locked;
assert(isfield(q,'name') && ~isempty(q.name),'MEG:Parameters','Parameter set needs a stable name.');
parameters=struct('name',char(q.name),'rho',q.rho,'tau',q.tau,'sigma',q.sigma);
end

function [manifest,diagnostics]=compute_cof(brainL,brainR,c,data,p,cfg,outputDir,branch)
valid=logical(data.reliable_frame_mask(:)');
[starts,ends]=logical_runs(valid);
parameters=struct('rho',p.rho,'tau',p.tau,'sigma',p.sigma);
chunks=struct('index',{},'filename',{},'interval_start',{},'interval_end',{},'bytes',{});
acc=init_diag(); chunkIndex=0;
for iSeg=1:numel(starts)
    frames=starts(iSeg):ends(iSeg);
    if numel(frames)<3,continue,end
    analyticL=analytic_mesh(brainL(:,:,frames),c.cii_mask_left);
    analyticR=analytic_mesh(brainR(:,:,frames),c.cii_mask_right);
    for first=1:cfg.cof_time_block_frames:numel(frames)-1
        last=min(first+cfg.cof_time_block_frames,numel(frames)); local=first:last;
        [uxL,uyL]=meg_cof_solve({real(analyticL(:,:,local)),imag(analyticL(:,:,local))}, ...
            c.cii_mask_left,parameters);
        [uxR,uyR]=meg_cof_solve({real(analyticR(:,:,local)),imag(analyticR(:,:,local))}, ...
            c.cii_mask_right,parameters);
        Ux=single(cat(2,uxL,uxR)); Uy=single(cat(2,uyL,uyR));
        input_frame_index=frames(local);
        interval_frame_index=input_frame_index(1:end-1);
        interval_time_seconds=data.time_seconds(interval_frame_index);
        chunkIndex=chunkIndex+1; filename=sprintf('velocity_chunk_%04d.mat',chunkIndex);
        save(fullfile(outputDir,filename),'Ux','Uy','input_frame_index', ...
            'interval_frame_index','interval_time_seconds','parameters','-v7.3');
        info=dir(fullfile(outputDir,filename));
        chunks(chunkIndex)=struct('index',chunkIndex,'filename',filename, ...
            'interval_start',interval_frame_index(1),'interval_end',interval_frame_index(end), ...
            'bytes',info.bytes);
        acc=update_diag(acc,Ux,Uy,analyticL(:,:,local),analyticR(:,:,local));
    end
end
if acc.n==0
    error('MEG:WorkflowNoValidCOFIntervals', ...
        'No contiguous reliable segment contained enough frames for COF.');
end
diagnostics=finish_diag(acc);
preprocessing=struct('mask',[c.cii_mask_left(:,1:end-1) c.cii_mask_right], ...
    'hemisphere_first',true,'hilbert_boundary','each contiguous reliable segment', ...
    'parameters',parameters,'maximum_imputed_hrf_mass',cfg.maximum_imputed_hrf_mass);
save(fullfile(outputDir,'cof_preprocessing.mat'),'preprocessing');
manifest=struct('schema_version',1,'status','complete','solver','meg_cof_solve.m', ...
    'branch',branch,'parameter_id',p.name,'parameters',parameters,'hemisphere_first',true, ...
    'frame_interval_seconds',cfg.output_tr_seconds,'input_frame_count',numel(valid), ...
    'reliable_input_frames',nnz(valid),'intervals_written',acc.n,'velocity_chunks',chunks, ...
    'reliable_frame_indices',find(valid), ...
    'diagnostics',diagnostics,'created',datestr(now,30));
write_json(fullfile(outputDir,'cof_manifest.json'),manifest);
end

function analytic=analytic_mesh(brain,mask)
[ny,nx,nt]=size(brain); active=find(mask(:)~=0);
flat=reshape(double(brain),ny*nx,nt); values=flat(active,:);
values=hilbert(values')';
out=complex(zeros(ny*nx,nt)); out(active,:)=values;
analytic=reshape(out,ny,nx,nt);
end

function a=init_diag()
a=struct('finite',0,'total',0,'speed',[],'rough_sum',0,'rough_n',0,'res_num',0,'res_den',0,'n',0);
end
function a=update_diag(a,ux,uy,imgL,imgR)
speed=sqrt(double(ux).^2+double(uy).^2); finite=isfinite(speed);
a.finite=a.finite+nnz(finite);a.total=a.total+numel(speed);
sample=sort(speed(finite)); if numel(sample)>5000,sample=sample(round(linspace(1,numel(sample),5000)));end
a.speed=[a.speed;sample(:)]; %#ok<AGROW>
dx=diff(double(ux),1,2);dy=diff(double(uy),1,1);
a.rough_sum=a.rough_sum+sum(abs(dx(isfinite(dx))))+sum(abs(dy(isfinite(dy))));
a.rough_n=a.rough_n+nnz(isfinite(dx))+nnz(isfinite(dy));
[rn,rd]=residual_ratio(imgL,imgR,ux,uy);a.res_num=a.res_num+rn;a.res_den=a.res_den+rd;
a.n=a.n+size(ux,3);
end
function [num,den]=residual_ratio(imgL,imgR,ux,uy)
num=0;den=0; images={imgL,imgR}; offset=0;
for h=1:2
    z=images{h}; M=size(z,1)-1;N=size(z,2)-1;
    dx=(z(1:M,2:N+1,1:end-1)-z(1:M,1:N,1:end-1)+z(2:M+1,2:N+1,1:end-1)-z(2:M+1,1:N,1:end-1))./2;
    dy=(z(2:M+1,1:N,1:end-1)-z(1:M,1:N,1:end-1)+z(2:M+1,2:N+1,1:end-1)-z(1:M,2:N+1,1:end-1))./2;
    dt=z(1:M,1:N,2:end)-z(1:M,1:N,1:end-1);
    cols=offset+(1:N); pred=dt+dx.*double(ux(:,cols,:))+dy.*double(uy(:,cols,:));
    f=isfinite(pred)&isfinite(dt);num=num+sum(abs(pred(f)).^2);den=den+sum(abs(dt(f)).^2);offset=offset+N;
end
end
function d=finish_diag(a)
d=struct('finite_fraction',a.finite/max(1,a.total),'speed_median',pct(a.speed,50), ...
    'speed_p99',pct(a.speed,99),'spatial_roughness',a.rough_sum/max(1,a.rough_n), ...
    'brightness_residual_ratio',a.res_num/max(eps,a.res_den),'intervals',a.n);
end
function v=pct(x,p)
if isempty(x),v=NaN;return,end;x=sort(x(:));pos=1+(numel(x)-1)*p/100;lo=floor(pos);hi=ceil(pos);v=x(lo)+(pos-lo)*(x(hi)-x(lo));
end
function row=make_row(subjectId,scanId,branch,p,d,outputDir)
row=struct('subject_id',subjectId,'scan_id',scanId,'branch',branch,'parameter_name',p.name, ...
    'rho',p.rho,'tau',p.tau,'sigma',p.sigma,'finite_fraction',d.finite_fraction, ...
    'speed_median',d.speed_median,'speed_p99',d.speed_p99,'spatial_roughness',d.spatial_roughness, ...
    'brightness_residual_ratio',d.brightness_residual_ratio,'intervals',d.intervals,'output_directory',outputDir);
end
function [s,e]=logical_runs(mask)
d=diff([false mask false]);s=find(d==1);e=find(d==-1)-1;
end
function branches=normalise_branches(value)
if ischar(value),branches={value};elseif isstring(value),branches=cellstr(value);else,branches=value;end
allowed={'exact_project1','bold_observed'};
if isempty(branches) || any(~ismember(branches,allowed)) || numel(unique(branches))~=numel(branches)
    error('MEG:WorkflowBranches','cof.branches must contain unique supported branch names.');
end
end
function validate_existing_manifest(manifest,branch,p,data)
required={'status','solver','branch','parameter_id','parameters','hemisphere_first', ...
    'input_frame_count','reliable_input_frames','reliable_frame_indices','diagnostics'};
for i=1:numel(required)
    if ~isfield(manifest,required{i}),error('MEG:WorkflowStaleManifest','Existing manifest lacks %s.',required{i});end
end
valid=logical(data.reliable_frame_mask(:)');
sameParameters=strcmp(char(manifest.parameter_id),p.name) && ...
    numeric_equal(manifest.parameters.rho,p.rho) && numeric_equal(manifest.parameters.tau,p.tau) && ...
    numeric_equal(manifest.parameters.sigma,p.sigma);
sameInput=double(manifest.input_frame_count)==numel(valid) && ...
    double(manifest.reliable_input_frames)==nnz(valid) && ...
    isequal(double(manifest.reliable_frame_indices(:))',double(find(valid)));
if ~strcmp(char(manifest.status),'complete') || ...
        ~strcmp(char(manifest.solver),'meg_cof_solve.m') || ...
        ~strcmp(char(manifest.branch),branch) || ~logical(manifest.hemisphere_first) || ...
        ~sameParameters || ~sameInput
    error('MEG:WorkflowManifestMismatch', ...
        'Existing COF manifest does not match the requested branch, parameters, or input support.');
end
end
function tf=numeric_equal(a,b)
tf=isfinite(double(a)) && isfinite(double(b)) && ...
    abs(double(a)-double(b))<=max(1e-15,16*eps(max(abs([double(a) double(b)]))));
end
function resultDir=cof_result_dir(cfg,parameterId,subjectId,scanId)
resultDir=fullfile(cfg.results_root,char(cfg.results_subdir),'cof',parameterId,subjectId,['scan-' scanId]);
end
function validate_metric_manifest(metricDir,p,cfg)
manifest=jsondecode(fileread(fullfile(metricDir,'divcurl_manifest.json')));
required={'status','method','parameters','velocity_parameters','chunks'};
for i=1:numel(required)
    if ~isfield(manifest,required{i}),error('MEG:WorkflowStaleMetrics','Metric manifest lacks %s.',required{i});end
end
q=manifest.velocity_parameters;s=manifest.parameters;
valid=strcmp(char(manifest.status),'complete') && strcmp(char(manifest.method),'meg_velocity_derivatives.m') && ...
    numeric_equal(q.rho,p.rho) && numeric_equal(q.tau,p.tau) && numeric_equal(q.sigma,p.sigma) && ...
    numeric_equal(s.filtSigma,cfg.metrics.filt_sigma) && ...
    logical(s.do_zscore)==logical(cfg.metrics.do_zscore) && ...
    logical(s.do_binarize)==logical(cfg.metrics.do_binarize) && ...
    numeric_equal(s.n_std,cfg.metrics.n_std) && numeric_equal(s.size_threshold,cfg.metrics.size_threshold);
if ~valid,error('MEG:WorkflowMetricMismatch','Existing metric manifest does not match the locked contract.');end
end
function write_json(filename,value)
fid=fopen(filename,'w');if fid<0,error('MEG:WorkflowJSONWrite','Cannot open %s.',filename);end
c=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(value),'char');clear c
end
