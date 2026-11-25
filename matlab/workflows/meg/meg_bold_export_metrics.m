function tableOut = meg_bold_export_metrics(configFile,subjectId,scanId,branch,parameterId,cofDir,metricDir)
%MEG_BOLD_EXPORT_METRICS Merge metric chunks and aggregate official blocks.

cfg=meg_bold_config(configFile);
manifest=jsondecode(fileread(fullfile(metricDir,'divcurl_manifest.json')));
chunks=manifest.chunks; tables=cell(numel(chunks),1);
for i=1:numel(chunks)
    s=load(fullfile(metricDir,chunks(i).filename),'PQ','PS','PQ_curl','PS_curl', ...
        'interval_frame_index','interval_time_seconds');
    n=numel(s.PQ);
    tables{i}=table(repmat({subjectId},n,1),repmat({scanId},n,1),repmat({branch},n,1), ...
        double(s.interval_frame_index(:)),double(s.interval_time_seconds(:)), ...
        double(s.PQ(:)),double(s.PS(:)),double(s.PQ_curl(:)),double(s.PS_curl(:)), ...
        'VariableNames',{'subject_id','scan_id','branch','interval_frame_index', ...
        'interval_time_seconds','PQ','PS','PQ_curl','PS_curl'});
end
tableOut=vertcat(tables{:});
resultDir=fullfile(cfg.results_root,char(cfg.results_subdir),'cof',parameterId,subjectId,['scan-' scanId]);
if ~isfolder(resultDir),mkdir(resultDir);end
writetable(tableOut,fullfile(resultDir,['frame_metrics_' branch '.csv']));

% Reuse the participant-specific official event decode audited in meg_native_conditions.
blockFile=fullfile(cfg.results_root,char(cfg.native_condition_subdir),subjectId,['scan-' scanId],'block_metrics.csv');
if isfile(blockFile)
    blocks=readtable(blockFile);
    rows=[];
    for i=1:height(blocks)
        in=tableOut.interval_time_seconds>=blocks.window_start_seconds(i) & ...
            tableOut.interval_time_seconds+cfg.output_tr_seconds<=blocks.window_end_seconds(i);
        r=table({subjectId},{scanId},{branch},blocks.block_index(i), ...
            blocks.condition_code(i),cell_text(blocks.condition,i), ...
            blocks.window_start_seconds(i),blocks.window_end_seconds(i),nnz(in), ...
            mean_or_nan(tableOut.PQ(in)),mean_or_nan(tableOut.PS(in)), ...
            'VariableNames',{'subject_id','scan_id','branch','block_index','condition_code', ...
            'condition','window_start_seconds','window_end_seconds','n_intervals','PQ','PS'});
        if isempty(rows),rows=r;else,rows=[rows;r];end %#ok<AGROW>
    end
    writetable(rows,fullfile(resultDir,['block_metrics_' branch '.csv']));
    if logical(cfg.metrics.save_maps_for_cross_modal)
        export_condition_maps(metricDir,chunks,blocks,resultDir,branch);
    end
else
    error('MEG:WorkflowOfficialBlocksMissing', ...
        'Participant-specific meg_native_conditions block table not found: %s',blockFile);
end
end

function export_condition_maps(metricDir,chunks,blocks,resultDir,branch)
conditionCodes=[3 1 2];
conditionNames={'fixation','0-back','2-back'};
patternNames={'Source','Sink','CCVortex','CVortex'};
sumMaps=[];counts=zeros(1,numel(conditionCodes));V_mask=[];
statsFile=fullfile(metricDir,'divcurl_statistics.mat');
if isfile(statsFile),tmp=load(statsFile,'V_mask');V_mask=tmp.V_mask;end
for iChunk=1:numel(chunks)
    s=load(fullfile(metricDir,chunks(iChunk).filename),'DivMap','CurlMap','interval_time_seconds');
    if ~isfield(s,'DivMap') || ~isfield(s,'CurlMap')
        error('MEG:WorkflowMotifMapsMissing','Metric chunks were not saved with D/C maps.');
    end
    if isempty(sumMaps),sumMaps=zeros(size(s.DivMap,1),size(s.DivMap,2),numel(conditionCodes),4);end
    maps={s.DivMap>0,s.DivMap<0,s.CurlMap>0,s.CurlMap<0};
    for iCondition=1:numel(conditionCodes)
        frameMask=false(numel(s.interval_time_seconds),1);
        conditionBlocks=blocks(blocks.condition_code==conditionCodes(iCondition),:);
        for b=1:height(conditionBlocks)
            frameMask=frameMask | (s.interval_time_seconds(:)>=conditionBlocks.window_start_seconds(b) & ...
                s.interval_time_seconds(:)<conditionBlocks.window_end_seconds(b));
        end
        counts(iCondition)=counts(iCondition)+nnz(frameMask);
        if any(frameMask)
            for iPattern=1:4
                sumMaps(:,:,iCondition,iPattern)=sumMaps(:,:,iCondition,iPattern)+ ...
                    sum(double(maps{iPattern}(:,:,frameMask)),3);
            end
        end
    end
end
condition_maps=nan(size(sumMaps),'single'); %#ok<NASGU>
for iCondition=1:numel(conditionCodes)
    if counts(iCondition)>0
        condition_maps(:,:,iCondition,:)=single(sumMaps(:,:,iCondition,:)./counts(iCondition));
    end
end
condition_frame_counts=counts; %#ok<NASGU>
Comparison={'2-back','0-back';'2-back','fixation';'0-back','fixation'}; %#ok<NASGU>
PatternNames=patternNames; %#ok<NASGU>
ConditionNames=conditionNames; %#ok<NASGU>
save(fullfile(resultDir,['condition_motif_maps_' branch '.mat']),'condition_maps', ...
    'condition_frame_counts','ConditionNames','PatternNames','Comparison','V_mask','-v7.3');
end
function v=mean_or_nan(x)
if isempty(x),v=NaN;else,v=mean(x(isfinite(x)));if isempty(v),v=NaN;end,end
end
function c=cell_text(column,i)
if iscell(column),c=column(i);elseif isstring(column),c={char(column(i))};else,c={char(string(column(i)))};end
end
