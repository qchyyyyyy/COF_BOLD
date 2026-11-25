function result = task_spatial_statistics(patterns,orderParameter,outputFile)
% Original latter-half spatial comparison block from visualize_tfMRI_HCP.m.
% patterns: cell(3,4), each grid x grid x participant for Fixation/Control/Task.
% orderParameter: cell(3,1), with the same axes.
assert(isequal(size(patterns),[3,4]) && numel(orderParameter)==3,'Unexpected map contract.');
pairs=[2 1;3 1;3 2];result=struct();
for k=1:3
    for m=1:4
        [~,p,~,stats]=ttest(patterns{pairs(k,1),m},patterns{pairs(k,2),m},'Dim',3);
        result.motif_p(:,:,m,k)=p;result.motif_t(:,:,m,k)=stats.tstat;
    end
    [~,p,~,stats]=ttest(orderParameter{pairs(k,1)},orderParameter{pairs(k,2)},'Dim',3);
    result.op_p(:,:,k)=p;result.op_t(:,:,k)=stats.tstat;
end
result.scope='Historical paired t-test formula on existing latter-half maps, without new correction';
if nargin>=3 && ~isempty(outputFile),assert(~isfile(outputFile),'Use a new output file.');save(outputFile,'result','-v7');end
end
