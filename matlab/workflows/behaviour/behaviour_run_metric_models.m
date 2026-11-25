function Model = behaviour_run_metric_models(configFile)
% Original extended-data size/order-parameter models, with explicit paths.
cfg = jsondecode(fileread(configFile));
assert(all(isfield(cfg,{'groupDir','behaviourCsv','outputDir'})), ...
    'Behaviour:Config','groupDir, behaviourCsv and outputDir are required.');
cfg.experiment='tfMRI_WM_LR'; cfg.session='WM'; cfg.method='COF';
cfg.hemispheres={'lh','rh'}; cfg.blockOrder=[8;1;6;4;5;7;2;3];
cfg.eventLength=404; cfg.blockColumns=1:8; cfg.useBoxCoxAccuracy=true;
cfg.conditionNames={'0-back','2-back'}; cfg.behaviors={'Acc','RT'};
cfg.behaviorLabels=containers.Map({'Acc','RT'},{'Box-Cox accuracy','Reaction time'});
cfg.appendixMetricDir=cfg.outputDir;
Data = behaviour_load_metrics(cfg);
Model = behaviour_metric_models(Data,cfg);
end
