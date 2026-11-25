function clinical_run(configFile)
%CLINICAL_RUN Export only the clinical summaries used in the manuscript.
cfg=jsondecode(fileread(configFile));
D=clinical_collect_cohort(cfg,cfg.dataset);
S=clinical_summarise(D);
clinical_output_directory(cfg.output_dir);
writetable(S.group_effects,fullfile(cfg.output_dir,'group_effects.csv'));
writetable(S.thresholds,fullfile(cfg.output_dir,'hc_thresholds.csv'));
writetable(S.proportions,fullfile(cfg.output_dir,'quadrant_proportions.csv'));
% This optional private plot-data file contains subject-level derived values.
if isfield(cfg,'save_plot_inputs') && cfg.save_plot_inputs
    save(fullfile(cfg.output_dir,'plot_inputs_private.mat'),'D','S','-v7.3');
end
if isfield(cfg,'disorganization_table') && ~isempty(cfg.disorganization_table)
    assert(strcmp(cfg.dataset,'ds000115'),'clinical:Dataset','Fig. 5g uses ds000115.');
    T=readtable(cfg.disorganization_table,'TextType','string');
    assert(numel(unique(T.participant_id))==height(T),'clinical:IDs','Duplicate historical table IDs.');
    [tf,idx]=ismember(D.subjectList,T.participant_id);
    assert(all(tf)&&height(T)==numel(tf),'clinical:IDs','Historical subject sets differ.');T=T(idx,:);
    expected=["HC-like efficient","Compensated low-PQ","Ineffective propagation","Under-propagation failure"];
    [tf,q]=ismember(T.Regime,expected);assert(all(tf),'clinical:Labels','Unknown frozen historical quadrant label.');
    assert(isequaln(double(T.Group01),D.group01) && isequaln([double(T.Age),double(T.Edu),double(T.Gender_01)],D.covMat),'clinical:Alignment','Historical covariates/group disagree.');
    trend=clinical_disorganization_trend(double(T.z_dis_4grp),q,D.group01,D.covMat);
    writetable(trend,fullfile(cfg.output_dir,'disorganization_trend.csv'));
end
copyfile(configFile,fullfile(cfg.output_dir,'settings.json'));
end
