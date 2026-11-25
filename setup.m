function root = setup(domain)
%SETUP Add the shared core and one scientific workflow to the MATLAB path.
% Run from the repository: setup('clinical'); clinical_run(configFile)
% Switching workflows removes only paths previously added by this function.
root=fileparts(mfilename('fullpath'));
persistent previousWorkflow
if ~isempty(previousWorkflow)
    old=strsplit(previousWorkflow,pathsep);
    current=strsplit(path,pathsep);
    for k=1:numel(old)
        if any(strcmp(old{k},current)),rmpath(old{k});end
    end
end
previousWorkflow='';
addpath(genpath(fullfile(root,'matlab','core')));
for name={'cifti_matlab','nanconv','raincloudplots','spm'}
    folder=fullfile(root,'third_party',name{1});
    if isfolder(folder),addpath(folder);end
end
if nargin>0 && strlength(string(domain))>0
    domain=validatestring(domain,{'simulations','rest','fmri','meg','behaviour','clinical','rnn','extended_data'});
    folder=fullfile(root,'matlab','workflows',domain);
    assert(isfolder(folder),'cof:Workflow','Workflow folder is missing.');
    previousWorkflow=genpath(folder);
    addpath(previousWorkflow);
end
end
