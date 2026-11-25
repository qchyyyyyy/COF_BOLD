function root=repository_root()
%REPOSITORY_ROOT Resolve this checkout independently of the working directory.
root=fileparts(mfilename('fullpath'));
while ~isfile(fullfile(root,'setup.m'))
    parent=fileparts(root);
    assert(~strcmp(parent,root) && ~isempty(parent),'cof:RepositoryRoot','Cannot locate repository setup.m.');
    root=parent;
end
end
