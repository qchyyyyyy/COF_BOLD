function simulate_wavefields(configFile)
% Frozen-input simulation export. Historical formulas are retained, not endorsed.
cfg=jsondecode(fileread(configFile));
assert(all(isfield(cfg,{'kind','seed','outputFile'})),'Explicit kind, seed and outputFile required');
assert(~isfile(cfg.outputFile),'Use a new output file');
rng(cfg.seed,'twister');payload=struct('kind',cfg.kind,'seed',cfg.seed,'panels',{{}});
switch cfg.kind
case 'packet'
    % Original exFig1_packet_propagation: rng(3) reset for every wave-count row.
    assert(cfg.seed==3,'Historical packet source fixes seed 3');
    [X,Y]=meshgrid(linspace(-6,6,100));t=linspace(1,2*pi*10,100);frames=1:8;
    for count=[1 2 5 10]
        rng(3,'twister');wave=zeros(100,100,100);K=randi(100,count,1)/100*.2;
        K(:,2)=sqrt(.2^2-K.^2);K=randi(10,count,1).*K;
        alpha=.1*(rand(size(wave))-.5);alpha=smooth3(alpha,'gaussian',5,10);
        for j=1:count
            kx=(1+1i*alpha).*K(j,1);ky=(1+1i*alpha).*K(j,2);
            amplitude=rand(1);phase=rand(1)*2*pi;
            wave=wave+amplitude.*exp(1i*(kx.*X+ky.*Y-reshape(t,1,1,100)+phase));
        end
        values=real(wave);limits=mean(values(:))+[-3 3]*std(values(:));
        payload.panels{end+1}=struct('name',sprintf('%d waves',count),'signal',values(:,:,frames),'K',K,'frames',frames,'limits',limits);
    end
case 'radial'
    % Historical analytic formulas, including their amplitude-gradient bias.
    [yi,xi]=meshgrid(-10:10);r=sqrt(xi.^2+yi.^2);T=100;c=1;omega=.5;
    data=zeros(21,21,T);
    for t=1:T,data(:,:,t)=cos(omega*(t-r/c))./(4*pi*r);end
    data(isinf(data))=nan;imaginary=reshape(imag(hilbert(reshape(data,[],T).')).',21,21,T);
    frames=2:7;v=nan(21,21,T);
    for t=1:T,v(:,:,t)=omega*imaginary(:,:,t)./(omega/c*imaginary(:,:,t)-data(:,:,t)./r);end
    radial={v,repmat(omega/c*ones(21),1,1,T),repmat((c*r.^2)./(1/omega^2+r.^2),1,1,T),repmat(c*ones(21),1,1,T)};
    names={'OF','PG (legacy descriptor)','COF (legacy formula)','Prescribed phase speed'};
    finite=data(isfinite(data));limits=mean(finite)+[-3 3]*std(finite);
    for j=1:4
        ux=xi./r.*radial{j};uy=yi./r.*radial{j};
        payload.panels{j}=struct('name',names{j},'signal',data(:,:,frames),'Ux',ux(:,:,frames),'Uy',uy(:,:,frames),'frames',frames,'limits',limits);
    end
case {'standing','gaussian'}
    T=100;times=linspace(1,2*pi*20,400);
    if strcmp(cfg.kind,'standing')
        [X,Y]=meshgrid(linspace(0,1,100));data=zeros(100,100,T);
        for t=1:T
            data(:,:,t)=sin(pi*X).*sin(2*pi*Y).*cos(.2*pi*sqrt(5)*times(t))+ ...
                sin(2*pi*X).*sin(pi*Y).*cos(.2*pi*sqrt(5)*times(t));
        end
        frames=26:4:46;tau=1e-3;
    else
        [X,Y]=meshgrid(linspace(-6,6,50),linspace(-6,6,100));data=zeros(100,50,T);
        envelope1=exp(-(X.^2+(Y-1.5).^2)/8);envelope2=exp(-(X.^2+(Y+1.5).^2)/8);
        decay=.5*cos(linspace(1,2*pi*2,400))+.5;
        for t=1:T
            data(:,:,t)=2*exp(-decay(t))*(cos(.5*times(t))*envelope1-sin(.5*times(t))*envelope2)+1e-5*randn(100,50);
        end
        frames=26:4:38;tau=1e-4;
        payload.seed_boundary='Original Gaussian script did not save a seed; configured seed is a reproducible reconstruction, not recovered original randomness.';
    end
    [M,N,~]=size(data);mask=ones(M,N);para=struct('tau',tau,'rho',1e-3,'sigma',1e-3);
    analytic=reshape(hilbert(reshape(data,[],T).').',M,N,T);
    for j=1:3
        if j==1,[ux,uy]=of_solve(data,mask,para);name='OF';
        elseif j==2,[ux,uy]=phase_gradient(angle(analytic));name='PG';
        else,[ux,uy]=cof_solve({data,imag(analytic)},mask,para);name='COF';end
        d=zeros(size(ux,1),size(ux,2),numel(frames));kernel=fspecial('gaussian',9,3);
        for k=1:numel(frames),d(:,:,k)=nanconv(divergence(ux(:,:,frames(k)),uy(:,:,frames(k))),kernel);end
        payload.panels{j}=struct('name',name,'signal',data(:,:,frames),'Ux',ux(:,:,frames),'Uy',uy(:,:,frames),'divergence',d,'frames',frames,'limits',[-2 2]);
    end
otherwise,error('Unknown simulation kind');
end
payload.job_id=getenv('SLURM_JOB_ID');payload.config=cfg;
parent=fileparts(cfg.outputFile);if ~isfolder(parent),mkdir(parent);end
save(cfg.outputFile,'payload','-v7');
end
