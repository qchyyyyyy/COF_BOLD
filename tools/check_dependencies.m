function report = check_dependencies()
%CHECK_DEPENDENCIES Exercise bundled IO and report optional external tools.
% Uses only synthetic temporary data. No participant data are read or written.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; oldRng = rng;
restoreState = onCleanup(@() restore_state(oldPath, oldRng)); %#ok<NASGU>
addpath(fullfile(root,'third_party','cifti_matlab'));
addpath(fullfile(root,'third_party','nanconv'));
addpath(fullfile(root,'third_party','raincloudplots'));
required = {'ft_read_cifti','ft_write_cifti','gifti','xmltree','nanconv'};
for k = 1:numel(required)
    resolved = which(required{k});
    assert(startsWith(resolved,fullfile(root,'third_party'),'IgnoreCase',ispc), ...
        'CPM:DependencyShadowing','%s resolves outside the bundled dependencies: %s',required{k},resolved);
end
report = struct('matlab',version,'bundled_io','not run','nanconv','not run', ...
    'raincloudplots','not run','spm_filter',which('spm_filter'), ...
    'spm_hrf',which('spm_hrf'),'fieldtrip',which('ft_defaults'), ...
    'boxcox',which('boxcox'));
if isempty(report.boxcox)
    report.financial_toolbox = 'not installed/on path: boxcox required for historical behaviour transformations';
else
    report.financial_toolbox = 'boxcox available (reported path); no substitute implementation used';
end
tmp = tempname; mkdir(tmp);
cleanupTmp = onCleanup(@() rmdir(tmp,'s')); %#ok<NASGU>

% Unequal numbers of vertices/time points detect accidental axis transposes.
s = struct();
s.pos = [(1:6)' zeros(6,2)];
s.brainstructure = [1;1;1;2;2;2];
s.brainstructurelabel = {'CORTEX_LEFT','CORTEX_RIGHT'};
s.time = (0:3)*0.72;
s.data = reshape(single(1:24),6,4)/8;
s.datadimord = 'pos_time';
filename = fullfile(tmp,'synthetic.dtseries.nii');
% The archived writer appends its type-specific suffix; pass an extensionless base.
ft_write_cifti(fullfile(tmp,'synthetic'),s,'parameter','data','writesurface',false);
r = ft_read_cifti(filename,'readdata',true,'readsurface',false);
assert(isequal(double(r.dtseries),double(s.data)), 'CPM:CiftiValues','CIFTI values/order changed.');
assert(isequal(r.brainstructure,s.brainstructure), 'CPM:CiftiStructure','CIFTI hemisphere order changed.');
assert(max(abs(r.time-s.time)) < 1e-7, 'CPM:CiftiTime','CIFTI sampling times changed.');

% A GZipBase64Binary surface exercises base64, XML and JVM inflate/deflate.
assert(usejava('jvm'),'CPM:JavaRequired','Bundled compressed GIfTI IO requires the MATLAB JVM.');
surface = struct('vertices',single([0 0 0;1 0 0;0 1 0]),'faces',int32([1 2 3]));
g = gifti(surface); surfaceFile = fullfile(tmp,'synthetic.surf.gii');
save(g,surfaceFile,'GZipBase64Binary');
g2 = gifti(surfaceFile);
assert(isequal(g2.vertices,g.vertices) && isequal(g2.faces,g.faces), ...
    'CPM:GiftiRoundTrip','Compressed GIfTI surface values/order changed.');
withSurfaces = ft_read_cifti(filename,'cortexleft',surfaceFile,'cortexright',surfaceFile);
assert(isequal(withSurfaces.dtseries,r.dtseries) && size(withSurfaces.pos,1)==6 ...
    && size(withSurfaces.tri,1)==2,'CPM:CiftiSurface','CIFTI-to-GIfTI surface loading failed.');
cpm_check_private_helpers();
report.bundled_io = 'passed: CIFTI, compressed GIfTI, XML and private warning helpers';

observed = nanconv([1 NaN 3],ones(1,3)/3,'nanout','1d');
assert(isnan(observed(2)) && all(isfinite(observed([1 3]))), ...
    'CPM:Nanconv','nanconv missing-value behaviour failed.');
report.nanconv = 'passed';

if exist('ksdensity','file') == 2 && license('test','Statistics_Toolbox')
    rng(90210,'twister');
    fig = figure('Visible','off');
    cleanupFigure = onCleanup(@() close(fig)); %#ok<NASGU>
    h = raincloud_plot([-2 -1 -.5 0 .2 .7 1 1.8 2.2], 'box_on',1);
    assert(~isempty(h) && all(cellfun(@(x) all(isgraphics(x)),h)), ...
        'CPM:Raincloud','Raincloud graphics were not created.');
    clf(fig);
    data = {[-2 -1 0 1 2];[-1.5 -.4 .5 1.8 2.5]};
    h2 = rm_raincloud(data,[0.2 0.4 0.7]);
    assert(isfield(h2,'p') && isgraphics(h2.p{1}), 'CPM:Raincloud','Repeated-measures raincloud failed.');
    report.raincloudplots = 'passed (default ks density)';
else
    report.raincloudplots = 'skipped: Statistics and Machine Learning Toolbox required';
end
if isempty(report.spm_filter) || isempty(report.spm_hrf)
    report.spm = 'not installed/on path: required for task-fMRI preprocessing';
else
    k = spm_filter(struct('RT',0.72,'row',1:80,'HParam',128));
    y = spm_filter(k,ones(80,1));
    hrf = spm_hrf(0.72);
    assert(all(isfinite(y)) && all(isfinite(hrf)) && abs(sum(hrf)-1)<1e-10, ...
        'CPM:SPM','External SPM filter/HRF smoke check failed.');
    report.spm = 'available; filter/HRF smoke check passed';
end
disp(report);
end

function restore_state(oldPath,oldRng)
path(oldPath); rng(oldRng);
end
