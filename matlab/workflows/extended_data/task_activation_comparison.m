function result = task_activation_comparison(time_r,seed,outputFile)
% Source: statistic_patternNum_corr.m, historical comparison block.
% time_r axes: four motifs x three methods (OF,PG,COF) x participants.
% The recovered source comments out the producer of time_r. Supply the original
% frozen array; this function does not claim that its upstream lineage is known.
assert(size(time_r,1)==4 && size(time_r,2)==3,'Expected 4 x 3 x participants.');
assert(all(isfinite(time_r(:))) && all(abs(time_r(:))<1),'Fisher transform requires finite abs(r)<1.');
rng(seed,'twister');comparisons=[1 2;2 3;1 3];p_values=zeros(3,4);
for i_pattern=1:4
    plot_data=abs(squeeze(time_r(i_pattern,:,:))');
    for k=1:size(comparisons,1)
        g1=comparisons(k,1);g2=comparisons(k,2);
        p_values(k,i_pattern)=permutation_mean_test(atanh(plot_data(:,g1)),atanh(plot_data(:,g2)),10000);
    end
end
result=struct('p_values',p_values,'comparisons',comparisons,'seed',seed, ...
    'scope','Original unadjusted comparison formula; historical seed and time_r provenance unresolved');
if nargin>=3 && ~isempty(outputFile),assert(~isfile(outputFile),'Choose a new output file.');save(outputFile,'result');end
end
