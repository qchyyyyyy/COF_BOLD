clear;
addpath(genpath('functions'));

% delete(gcp('nocreate'));
% partool = parpool('local',35);
% partool.IdleTimeout = Inf;

OutputPath = 'data/results/';
PE = {'LR','RL'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = {'lh','rh'};

for pe = 1
    if pe==1
        Condition = {'MOTOR','LANGUAGE','RELATIONAL','WM','REST1'};% 'WM','RELATIONAL',
    else
        Condition = {'RELATIONAL','WM'};
    end

    for cdn = 1:length(Condition)
        if strcmp(Condition{cdn}(1:end-1),'REST')
            DataPath = ['data/fMRI/rfMRI_REST1_' PE{pe} '/'];
        else
            DataPath = ['data/fMRI/tfMRI_' Condition{cdn} '_' PE{pe} '/'];
        end
        
        
        File = dir(fullfile(DataPath));
        SubjectList = cell(size(File));
        for m = 1:length(File)
            SubjectList{m}=File(m).name;
        end
        SubjectList(1:2)=[];
        save(['data/results/' Condition{cdn} '_' PE{pe} '_SubList.mat'],'SubjectList');
        for hs = 1:2
            switch hemisphere{hs}
                case 'lh'
                    [mask, label, oi_x, oi_y, V_mask, mask1] = masklabel(parcellation,'L');
                case 'rh'
                    [mask, label, oi_x, oi_y, V_mask, mask1] = masklabel(parcellation,'R');
            end
            
            NumSource = [];NumSink = [];NumCCVortex = [];NumCVortex = [];
            SourceMap = [];
            SinkMap = [];
            CCVortexMap = [];
            CVortexMap = [];
            
            for sub = length(SubjectList):-1:1
                disp(['sub = ' num2str(sub)]);
                load([DataPath,SubjectList{sub}, '/VelField_' hemisphere{hs} '_hilbert.mat']);
                
                Ux = VelField.Ux;
                Uy = VelField.Uy;
                T = size(Ux,3);
                
                filtSigma = 3;
                filtWidth = round(3*filtSigma);
                imageFilter=fspecial('gaussian',filtWidth,filtSigma);
                DivMap=zeros(size(Ux));CurlMap = DivMap;
                for t = T:-1:1
                    U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);

                    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
                    div_oneframe(V_mask==0)=nan;
                    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
                    
                    [curl_oneframe,~]= curl(Ux(:,:,t)./U,Uy(:,:,t)./U);
                    curl_oneframe(V_mask==0)=nan;
                    CurlMap(:,:,t) = nanconv(curl_oneframe,imageFilter);
                end
                
                CurlMap(repmat(V_mask,1,1,T)==0)=nan;
                DivMap(repmat(V_mask,1,1,T)==0)=nan;
                
                k=2;
                DivMean = mean(DivMap,'all','omitnan');
                DivStd = std(DivMap,[],'all','omitnan');
                DivMap(DivMap<DivMean+k*DivStd & DivMap>DivMean-k*DivStd)=0;
                DivMap(DivMap>0)=1;
                DivMap(DivMap<0)=-1;
                
                CurlMean = mean(CurlMap,'all','omitnan');
                CurlStd = std(CurlMap,[],'all','omitnan');
                CurlMap(CurlMap<CurlMean+k*CurlStd & CurlMap>CurlMean-k*CurlStd)=0;
                CurlMap(CurlMap>0)=1;
                CurlMap(CurlMap<0)=-1;
                
                SoM = DivMap;SoM(SoM<0)=0;
                SiM = -DivMap;SiM(SiM<0)=0;
                CCM = CurlMap;CCM(CCM<0)=0;
                CM = -CurlMap;CM(CM<0)=0;
                
                SourceMap(:,:,sub) = mean(SoM,3,'omitnan');
                SinkMap(:,:,sub) = mean(SiM,3,'omitnan');
                CCVortexMap(:,:,sub) = mean(CCM,3,'omitnan');
                CVortexMap(:,:,sub) = mean(CM,3,'omitnan');
                
                
                parfor t = 1:T
                    [NumSource(sub,t), NumSink(sub,t)] = patternNum(DivMap(:,:,t));
                    [NumCCVortex(sub,t), NumCVortex(sub,t)] = patternNum(CurlMap(:,:,t));
                end
                
            end

            
            save([OutputPath 'numPattern_results/voxel/' hemisphere{hs} '_numPattern_individual_' Condition{cdn} '_' PE{pe} '_voxel.mat'],...
                'SourceMap','SinkMap',...
                'CCVortexMap','CVortexMap','-v7.3');
            
            save([OutputPath 'numPattern_results/global/' hemisphere{hs} '_numPattern_individual_' Condition{cdn} '_' PE{pe} '.mat'],...
                'NumSource','NumSink',...
                'NumCCVortex','NumCVortex','-v7.3');
        end
    end
    
end
%% sub function
function [numpattern1, numpattern2] = patternNum(Map)
ind = find(Map==1);
if length(ind)<10
    numpattern1=0;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern1 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern1);
    ps_k=zeros(numpattern1,1);
    for k = numpattern1:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    numpattern1 = numpattern1 - length(find(ps_k<10));
end

ind = find(Map==-1);
if length(ind)<10
    numpattern2=0;
else
    [I,J] = ind2sub(size(Map),ind);
    Z = linkage([I J],'single');
    numpattern2 = length(find(Z(:,3)>sqrt(2)))+1;
    Cluster_ind = cluster(Z,'MaxClust',numpattern2);
    ps_k=zeros(numpattern2,1);
    for k = numpattern2:-1:1
        ps_k(k) = length(find(Cluster_ind==k));
    end
    numpattern2 = numpattern2 - length(find(ps_k<10));
end

end
