function T = behaviour_export_inputs(configFile)
% Export the exact participant-level inputs for the manuscript global models.
% Accuracy is Box-Cox transformed across all selected blocks before averaging.
cfg = jsondecode(fileread(configFile));
assert(all(isfield(cfg,{'group_dir','behaviour_csv','output_file'})), ...
    'Behaviour:Config','group_dir, behaviour_csv and output_file are required.');
assert(~isfile(cfg.output_file),'Behaviour:OutputExists','Use a new output file.');
if isfield(cfg,'historical_input_hashes') && ~isempty(cfg.historical_input_hashes)
    behaviour_verify_inputs(cfg.group_dir,cfg.historical_input_hashes);
end
p = 'tfMRI_WM_LR';
subjects = behaviour_subject_list(fullfile(cfg.group_dir,[p '_SubList.mat']));
I = load(fullfile(cfg.group_dir,[p '_info.mat']),'Acc','RT','Stim_Type_subs');
ids = str2double(string(subjects(:)));
assert(all(isfinite(ids)) && numel(unique(ids))==numel(ids), ...
    'Behaviour:Subjects','Participant IDs must be finite and unique.');
Acc = I.Acc;
Acc(Acc==0) = nan;
sel = sum((squeeze(I.Stim_Type_subs(:,3,:))-[8;1;6;4;5;7;2;3]).^2)==0 & ~isnan(sum(Acc,2))';
ix = find(sel(:));
assert(~isempty(ix),'Behaviour:Subjects','No participants meet the original inclusion rule.');
if isfield(cfg,'expected_n')
    assert(numel(ix)==cfg.expected_n,'Behaviour:Subjects','Unexpected selected cohort size.');
end
if isfield(cfg,'historical_subject_csv') && ~isempty(cfg.historical_subject_csv)
    H = readtable(cfg.historical_subject_csv);
    assert(isequal(ids(ix),H.Subject(:)),'Behaviour:Order','Historical participant order differs.');
end
timing = mean(I.Stim_Type_subs(:,:,sel),3);
timing(:,2) = sum(timing(:,1:2),2);
ev = zeros(404,1);
for k = 1:8
    frame = round(timing(k,1)):round(timing(k,2))-1;
    frame = frame(frame>=1 & frame<=404);
    frame = frame(1:min(38,numel(frame)));
    ev(frame) = timing(k,3);
end
assert(all(arrayfun(@(k)nnz(ev==k)>0,1:8)), ...
    'Behaviour:Timing','Each of the eight task blocks needs valid frames.');
A = Acc(ix,1:8);
R = I.RT(ix,1:8);
L = load(fullfile(cfg.group_dir,[p '_numPattern_global_lh_COF.mat']),'NumSource','NumSink');
U = load(fullfile(cfg.group_dir,[p '_numPattern_global_rh_COF.mat']),'NumSource','NumSink');
assert(isequal(size(L.NumSource),size(U.NumSource),size(L.NumSink),size(U.NumSink)), ...
    'Behaviour:Shape','Bilateral source/sink array sizes differ.');
assert(size(L.NumSource,1)==numel(ids) && size(L.NumSource,2)==404, ...
    'Behaviour:Shape','Expected one complete 404-frame row per participant.');
assert(exist('boxcox','file')==2,'Behaviour:Dependency', ...
    'The historical Box-Cox calculation requires the original MATLAB boxcox function.');
[bc,lambda] = boxcox(A(:));
bc = reshape(bc,size(A));
N = L.NumSource(ix,:)+U.NumSource(ix,:);
K = L.NumSink(ix,:)+U.NumSink(ix,:);
B = readtable(cfg.behaviour_csv);
assert(numel(unique(B.Subject))==height(B),'Behaviour:Subjects','Behaviour IDs are duplicated.');
[ok,bi] = ismember(ids(ix),B.Subject);
assert(all(ok),'Behaviour:Alignment','Behaviour CSV lacks a selected participant.');
B = B(bi,:);
age = nan(height(B),1);
for k = 1:height(B)
    v = str2double(split(string(B.Age(k)),'-'));
    v = v(isfinite(v));
    % Retain the historical HCP open-ended-age fallback.
    if isempty(v), age(k)=36; else, age(k)=mean(v); end
end
T = table(ids(ix),age,double(string(B.Gender)=="M"), ...
    'VariableNames',{'Subject','age','sex'});
for c = [0 2]
    if c==0, frames=find(ev>0 & ev<5); blocks=1:4;
    else, frames=find(ev>4); blocks=5:8; end
    suffix = num2str(c);
    T.(['source_' suffix]) = mean(N(:,frames),2);
    T.(['sink_' suffix]) = mean(K(:,frames),2);
    T.(['accuracy_' suffix]) = mean(A(:,blocks),2);
    T.(['boxcox_' suffix]) = mean(bc(:,blocks),2);
    T.(['rt_' suffix]) = mean(R(:,blocks),2);
end
if exist('H','var')
    assert(max(abs(T.source_2-H.PQ_source_2back))<1e-10);
    assert(max(abs(T.sink_2-H.PQ_sink_2back))<1e-10);
    assert(max(abs(T.boxcox_2-H.Acc_2back_boxcox))<1e-10);
    assert(max(abs(T.accuracy_2-H.Acc_2back_raw))<1e-10);
end
parent = fileparts(cfg.output_file);
if ~isempty(parent) && ~isfolder(parent), mkdir(parent); end
writetable(T,cfg.output_file);
meta = struct('participants',height(T),'boxcox_lambda',lambda, ...
    'two_back_frames',nnz(ev>4),'zero_back_frames',nnz(ev>0 & ev<5));
fid = fopen([cfg.output_file '.json'],'w');
assert(fid>=0,'Behaviour:Output','Cannot open provenance output.');
cleanup = onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(meta));
end
