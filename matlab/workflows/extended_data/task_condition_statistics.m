function result = task_condition_statistics(configFile)
% Original ED8/ED9 condition aggregation and independent permutation tests.
% cfg: group_dir, experiment, output_file, seed, branch ('ED8' or 'ED9').
c=jsondecode(fileread(configFile));rng(c.seed,'twister');
assert(~isfile(c.output_file),'Use a new output file.');
E=load(fullfile(c.group_dir,[c.experiment '_EVs.mat']),'Events','Stim_Type_subs');
Events=E.Events;keep=true(size(Events,1),1);
if strcmp(c.branch,'ED9')
    block_setting=squeeze(E.Stim_Type_subs(:,3,:));
    keep=(sum((block_setting-mean(block_setting,2)).^2)<1)';
    Events=repmat(mean(Events(keep,:),1),nnz(keep),1);patterns=1:4;
else
    assert(strcmp(c.branch,'ED8'),'branch must be ED8 or ED9');patterns=1:2;
end
names={'Source','Sink','CCVortex','CVortex'};metrics={'numPattern','sizePattern','orderParam'};
result=struct('branch',c.branch,'seed',c.seed,'retained_n',nnz(keep),'permutations',5000,'correction','none');
for i_metric=1:numel(metrics)
    key=metrics{i_metric};S=load(fullfile(c.group_dir,[c.experiment '_' key '_global_LR_COF.mat']));
    if strcmp(key,'orderParam'),TC=S.order_param_global(keep,:);
    else
        pre='Num';if strcmp(key,'sizePattern'),pre='Size';end
        TC=[];for k=patterns,TC(:,:,k)=S.([pre names{k}])(keep,:);end
    end
    Time_included=true(1,size(TC,2));Time_included([1:10 end-9:end])=false;
    data=cell(1,3);
    for i_sub=1:size(TC,1)
        for i_pattern=1:size(TC,3)
            tmp_TC=squeeze(TC(i_sub,:,i_pattern));
            data{1}(i_sub,i_pattern)=mean(tmp_TC(Events(i_sub,:)==0 & Time_included));
            data{2}(i_sub,i_pattern)=mean(tmp_TC(Events(i_sub,:)==1 & Time_included));
            data{3}(i_sub,i_pattern)=mean(tmp_TC(Events(i_sub,:)==2 & Time_included));
        end
    end
    p_values=[];t_values=[];
    for i_pattern=1:size(TC,3)
        [p_values(i_pattern,1),t_values(i_pattern,1)]=permutation_mean_test(data{2}(:,i_pattern),data{1}(:,i_pattern),5000);
        [p_values(i_pattern,2),t_values(i_pattern,2)]=permutation_mean_test(data{3}(:,i_pattern),data{1}(:,i_pattern),5000);
        [p_values(i_pattern,3),t_values(i_pattern,3)]=permutation_mean_test(data{3}(:,i_pattern),data{2}(:,i_pattern),5000);
    end
    result.(key)=struct('condition_values',{data},'p_values',p_values,'t_values',t_values);
end
result.scope='Historical independent permutation formula; no family adjustment; seed is explicit reconstruction';
save(c.output_file,'result','-v7');
end
