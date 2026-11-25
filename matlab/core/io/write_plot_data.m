function filename=write_plot_data(renderer,args,prefixIndex)
% Save computed plot inputs without invoking a renderer. Outputs may contain restricted
% derived values; store/transfer privately and never include them in Git.
assert(~isempty(regexp(renderer,'^[A-Za-z]\w*$','once')),'Invalid renderer');
prefix=char(args{prefixIndex}); filename=[prefix '.render.mat'];
assert(~isfile(filename),'cof:OutputExists','Choose a fresh plot-data output path.');
parent=fileparts(filename);if ~isfolder(parent),mkdir(parent);end
schema_version=1;
save(filename,'renderer','args','prefixIndex','schema_version','-v7.3');
end
