function surrogate_historical_summary(configFile)
% Prepared replay of recovered historical sources; requires original inputs.
% Cluster compute/export only. Do not run the archived mixed plotting scripts.
c=jsondecode(fileread(configFile));
out=c.output_dir;assert(~isfolder(out),'Use a new output directory.');mkdir(out);
C=load(c.coordinate_file,'cii_mask_left','cii_mask_right');
vl=conv2(double(C.cii_mask_left==1),ones(2)/4,'valid')==1;
vr=conv2(double(C.cii_mask_right==1),ones(2)/4,'valid')==1;mask=[vl vr];
if isfield(c,'hist_real_file') && ~isempty(c.hist_real_file)
    assert(isfile(c.hist_real_file)&&isfile(c.hist_surrogate_2D_file),'Original matched real/2D-surrogate input pair required');
    files={c.hist_real_file,c.hist_surrogate_2D_file};edges_div=-2:.05:2;edges_curl=-2:.01:2;
    prob_div=zeros(2,numel(edges_div)-1);prob_curl=zeros(2,numel(edges_curl)-1);counts=zeros(2,2);
    for k=1:2
        S=load(files{k},'data');x=S.data.BrainImg_lh(:,:,201:301);
        % Preserve actual positional arguments, not misleading caller names:
        % old call (tau,rho,sigma) binds to solver signature (rho,tau,sigma).
        [ux,uy]=surrogate_estimate_velocity(x,C.cii_mask_left,5e-4,1e-8,8);
        [d,w]=surrogate_velocity_derivatives(ux(:,:,11:90),uy(:,:,11:90),vl,3,0);
        valid=repmat(vl,1,1,80);a=d(valid);b=w(valid);
        prob_div(k,:)=histcounts(a,edges_div,'Normalization','probability');
        prob_curl(k,:)=histcounts(b,edges_curl,'Normalization','probability');
        counts(k,:)=[nnz(isfinite(a)),nnz(isfinite(b))];
        raw_div{k}=a;raw_curl{k}=b;
    end
    % Exact historical pooled vartest2 formula: descriptive source replay only.
    [~,p_div,ci_div,stats_div]=vartest2(raw_div{1},raw_div{2});
    [~,p_curl,ci_curl,stats_curl]=vartest2(raw_curl{1},raw_curl{2});
    save(fullfile(out,'surrogate_pooled_variance_tests.mat'),'p_div','ci_div','stats_div','p_curl','ci_curl','stats_curl');
    save(fullfile(out,'surrogate_histogram.mat'),'prob_div','prob_curl','edges_div','edges_curl','counts','-v7');
end
if isfield(c,'density_lh_file') && ~isempty(c.density_lh_file)
    files={c.density_lh_file,c.density_rh_file};names={'SourceMap','SinkMap','CCVortexMap','CVortexMap'};
    maps=zeros(200,500,4);ns=zeros(2,4);
    for h=1:2
        assert(isfile(files{h}),'Original hemisphere density file required');m=matfile(files{h});
        for k=1:4
            shape=size(m,names{k});assert(isequal(shape(1:2),[200 250]));ns(h,k)=shape(3);total=zeros(200,250);
            for first=1:32:shape(3)
                last=min(first+31,shape(3));block=m.(names{k})(:,:,first:last);total=total+sum(double(block),3);
            end
            maps(:,(h-1)*250+(1:250),k)=total/shape(3);
        end
    end
    assert(all(ns(:)==ns(1)),'Hemisphere/motif sample counts differ');n=ns(1);
    save(fullfile(out,'surrogate_density.mat'),'maps','mask','n','-v7');
end
end
