function clinical_aggregate_pq(configFile)
%CLINICAL_AGGREGATE_PQ Assemble participants in an explicit, checked row order.
% cfg.subject_ids (unique); cfg.run_dirs is participants x runs, each a folder
% produced by clinical_compute. ds000115 order: 0-, 1-, 2-back, 475 frames each.
% Normalize derivatives over all runs together BEFORE counting motifs.
cfg=jsondecode(fileread(configFile)); ids=string(cfg.subject_ids(:));
assert(numel(unique(ids))==numel(ids),'clinical:IDs','Duplicate subject IDs.');
if iscell(cfg.run_dirs) && all(cellfun(@iscell,cfg.run_dirs(:)))
    rows=cellfun(@(r)reshape(string(r),1,[]),cfg.run_dirs(:),'UniformOutput',false);
    assert(~isempty(rows) && numel(unique(cellfun(@numel,rows)))==1,'clinical:Shape','Run lists must have equal lengths.');
    files=vertcat(rows{:});
else
    files=string(cfg.run_dirs);
    if numel(ids)==1,files=reshape(files,1,[]);end
end
assert(size(files,1)==numel(ids)&&size(files,2)>0,'clinical:Shape','One nonempty run-list row is required per participant.');
m=load(fullfile(cfg.template_root,'coordinate_cii_v2.mat'),'cii_mask_left','cii_mask_right');
ml=m.cii_mask_left;mr=m.cii_mask_right;ml(isnan(ml))=0;mr(isnan(mr))=0;
ml=conv2(ml,ones(2)/4,'valid');mr=conv2(mr,ones(2)/4,'valid');
ml(ml~=1)=0;mr(mr~=1)=0;mask=cat(2,ml,mr);
method=validatestring(cfg.method,{'OF','PG','COF'});
NumSource=[];NumSink=[];
for i=1:numel(ids)
    left=cell(1,size(files,2));right=left;
    for j=1:size(files,2)
        l=load(fullfile(files(i,j),['VelField_lh_' method '.mat']),'VelField');
        r=load(fullfile(files(i,j),['VelField_rh_' method '.mat']),'VelField');
        left{j}=l.VelField;right{j}=r.VelField;
        if isfield(cfg,'frames_per_run')
            expected=cfg.frames_per_run+strcmp(method,'PG');
            assert(size(l.VelField.Ux,3)==expected&&size(r.VelField.Ux,3)==expected,'clinical:Frames','Unexpected run length.');
        end
    end
    x=clinical_count_pq(left,right,mask,method);src=x.NumSource;snk=x.NumSink;
    if i>1,assert(numel(src)==size(NumSource,2),'clinical:Frames','Participants have different run lengths.');end
    NumSource(i,:)=src;NumSink(i,:)=snk; %#ok<AGROW>
end
clinical_output_directory(cfg.output_dir);
subjectList=ids; %#ok<NASGU>
[~,name,ext]=fileparts(cfg.output_name);assert(strcmp(cfg.output_name,[name ext])&&strcmp(ext,'.mat'),'clinical:Output','output_name must be a MAT basename.');
save(fullfile(cfg.output_dir,cfg.output_name),'NumSource','NumSink','subjectList');
end
