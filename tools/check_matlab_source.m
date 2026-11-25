function report=check_matlab_source(outputFile)
%CHECK_MATLAB_SOURCE Run MATLAB Code Analyzer on author-owned source files.
root=fileparts(fileparts(mfilename('fullpath')));
files=[dir(fullfile(root,'matlab','**','*.m'));dir(fullfile(root,'tests','matlab','*.m'));dir(fullfile(root,'tools','*.m'));dir(fullfile(root,'setup.m'))];
issues=struct('file',{},'id',{},'message',{},'line',{});syntaxErrors=0;
for k=1:numel(files)
    path=fullfile(files(k).folder,files(k).name);messages=checkcode(path,'-id');
    for j=1:numel(messages)
        m=messages(j);issues(end+1)=struct('file',erase(path,[root filesep]),'id',m.id,'message',m.message,'line',m.line); %#ok<AGROW>
        if any(strcmp(m.id,{'PARSE','SYNER','ENDCT','ENDCT0','UNBAL','ERRROW'})),syntaxErrors=syntaxErrors+1;end
    end
end
report=struct('files',numel(files),'syntax_errors',syntaxErrors,'issues',issues);
if nargin>0,fid=fopen(outputFile,'w');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));end %#ok<NASGU>
assert(syntaxErrors==0,'cof:Syntax','MATLAB parse errors found.');
fprintf('Code Analyzer: %d author files, %d syntax errors; %d total diagnostics.\n',numel(files),syntaxErrors,numel(issues));
end
