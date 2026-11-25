function filename = meg_subject_manifest(cfg, projectRoot)
% Explicit external participant manifest; never bundle HCP identifiers.
filename=getenv('COF_MEG_SUBJECT_MANIFEST');
if isempty(filename)
    filename=fullfile(projectRoot,char(cfg.subject_manifest_relative));
end
assert(isfile(filename),'NC:MEGManifest','Set COF_MEG_SUBJECT_MANIFEST to the locked cohort list.');
end
