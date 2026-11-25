function clinical_compute(configFile)
%CLINICAL_COMPUTE Preprocess declared scans and estimate VVFs.
% cfg.runs is an explicit scan list; each has input_cifti, output_dir, tr,
% and kind ('task' or 'rest'). Do not infer cohorts from a directory listing.
cfg=jsondecode(fileread(configFile));
assert(all(isfield(cfg,{'runs','template_root','methods','rho','tau','sigma'})),'clinical:Config','Missing compute settings.');
old=getenv('COF_TEMPLATE_ROOT');restore=onCleanup(@()setenv('COF_TEMPLATE_ROOT',old)); %#ok<NASGU>
setenv('COF_TEMPLATE_ROOT',cfg.template_root);
m=load(fullfile(cfg.template_root,'coordinate_cii_v2.mat'),'cii_mask_left','cii_mask_right');
methods=cellstr(string(cfg.methods));
assert(all(ismember(string(methods),["OF","PG","COF"])),'clinical:Method','Unsupported estimator.');
for k=1:numel(cfg.runs)
    run=cfg.runs(k);assert(isfile(run.input_cifti),'clinical:Input','Missing scan: %s',run.input_cifti);
    clinical_output_directory(run.output_dir);
    switch run.kind
        case 'task',ok=preprocess_task(run.input_cifti,run.output_dir,run.tr);
        case 'rest',ok=preprocess_rest(run.input_cifti,run.output_dir,run.tr);
        otherwise,error('clinical:Kind','Use task or rest.');
    end
    assert(ok,'clinical:Preprocessing','Preprocessing failed: %s',run.input_cifti);
    for j=1:numel(methods)
        method=methods{j}; input=struct('path',run.output_dir,'fileName','TC_mesh_r_f_zscore.mat','hemisphere','lh');
        estimate_velocity(input,m.cii_mask_left,method,cfg.rho,cfg.tau,cfg.sigma);
        input.hemisphere='rh';estimate_velocity(input,m.cii_mask_right,method,cfg.rho,cfg.tau,cfg.sigma);
    end
    settings=run;save(fullfile(run.output_dir,'scan_settings.mat'),'settings','cfg');
end
end
