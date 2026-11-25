function path = template_path(relative)
% Resolve externally installed template assets without changing the working folder.
root=getenv('COF_TEMPLATE_ROOT');
if isempty(root)
    dataRoot=getenv('COF_DATA_ROOT');
    if ~isempty(dataRoot), root=fullfile(dataRoot,'templates'); end
end
assert(~isempty(root),'release:TemplateRoot','Set COF_TEMPLATE_ROOT to the external template directory.');
path=fullfile(root,relative);
end
