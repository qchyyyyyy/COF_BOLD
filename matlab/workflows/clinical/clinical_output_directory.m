function clinical_output_directory(folder)
%CLINICAL_OUTPUT_DIRECTORY Require a new external output directory.
f=java.io.File(folder);assert(f.isAbsolute(),'clinical:Output','Output must be an absolute path.');
dest=char(f.getCanonicalPath());root=char(java.io.File(repository_root()).getCanonicalPath());
assert(~strcmpi(dest,root)&&~startsWith(lower(dest),[lower(root) filesep]),'clinical:Output','Store participant-derived outputs outside the checkout.');
assert(~isfolder(dest) && ~isfile(dest),'clinical:OutputExists','Choose a new output directory: %s',dest);
[ok,msg]=mkdir(dest);assert(ok,'clinical:Output','Cannot create output directory: %s',msg);
end
