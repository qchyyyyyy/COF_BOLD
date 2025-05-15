clear;
addpath(genpath('functions'));

Condition = {'WM','RELATIONAL'};
PE = {'LR'};
parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';

[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere(1));
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere(2));

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
%%
for cdn = 1
    for pe = 1
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
        
        
        load([Condition{cdn} '_' PE{pe} '_info_new.mat']);
        switch Condition{cdn}
            case 'WM'
                T = 404;
            case 'RELATIONAL'
                T = 231;
        end
        SubjectList(BadInd)=[];
        
        Events = zeros(T,1);
        Time = 10:T-9;
        
        for bb = 1:8
            idx = ceil(meanStimType(bb,1)-0.5):fix(meanStimType(bb,2)-0.5);
            idx(idx>T)=[];idx(idx<1)=[];
            Events(idx) = meanStimType(bb,3);
        end
        
        Events(Events>0 & Events<5)=1;
        Events(Events>4)=2;
        Events = Events+1;
        
        
        SourceMap=zeros([size(V_mask),T]);
        SinkMap=SourceMap;
        CCVortexMap=SourceMap;
        CVortexMap=SourceMap;
        
        for sub = length(SubjectList):-1:1
            disp(['sub = ' num2str(sub)]);
            L = load([DataPath,SubjectList{sub}, '/VelField_lh_hilbert.mat']);
            R = load([DataPath,SubjectList{sub}, '/VelField_rh_hilbert.mat']);
            
            Ux = cat(2,L.VelField.Ux,R.VelField.Ux);
            Uy = cat(2,L.VelField.Uy,R.VelField.Uy);clearvars L R
            %             T = size(Ux,3);
            
            filtSigma = 3;
            filtWidth = round(3*filtSigma);
            imageFilter=fspecial('gaussian',filtWidth,filtSigma);
            
            DivMap=[];CurlMap=[];
            for t = T:-1:1
                U = sqrt(Ux(:,:,t).^2+Uy(:,:,t).^2);
                %         div_oneframe= divergence(Ux(:,:,t)./U,Uy(:,:,t)./U);
                div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
                div_oneframe(V_mask==0)=nan;
                DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
                
                %             [curl_oneframe,~]= curl(Ux(:,:,t),Uy(:,:,t));
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

            Feature = [map2vec(SoM(:,:,Time),V_mask);...
                map2vec(SiM(:,:,Time),V_mask);...
                map2vec(CCM(:,:,Time),V_mask);...
                map2vec(CM(:,:,Time),V_mask)]';
            Feature = zscore(Feature);Cov = corr(Feature','rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_all(:,sub) = calc_decoding_acc(idx,K,Events,Time);
            
            Feature = [map2vec(SoM(:,:,Time),V_mask);...
                map2vec(SiM(:,:,Time),V_mask)]';
            Feature = zscore(Feature);Cov = corr(Feature','rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_div(:,sub) = calc_decoding_acc(idx,K,Events,Time);
            
            Feature = [map2vec(CCM(:,:,Time),V_mask);...
                map2vec(CM(:,:,Time),V_mask)]';
            Feature = zscore(Feature);Cov = corr(Feature','rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_curl(:,sub) = calc_decoding_acc(idx,K,Events,Time);
            
            SourceMap = SourceMap + SoM;
            SinkMap = SinkMap + SiM;
            CCVortexMap = CCVortexMap + CCM;
            CVortexMap = CVortexMap + CM;
        end
        
        SourceMap = SourceMap./length(SubjectList);
        SinkMap = SinkMap./length(SubjectList);
        CCVortexMap = CCVortexMap./length(SubjectList);
        CVortexMap = CVortexMap./length(SubjectList);
        
        save(['data/results/numPattern_results/trial/' hemisphere '_numPattern_individual_' Condition{cdn} '_' PE{pe} '_trial.mat'],...
            'SourceMap','SinkMap',...
            'CCVortexMap','CVortexMap','-v7.3');
        save(['data/results/decoding_results/' hemisphere '_decodingAcc_Pattern_individual_' Condition{cdn} '_' PE{pe} '.mat'],...
            'decodingAcc_all','decodingAcc_div',...
            'decodingAcc_curl','-v7.3');
    end
end
%% sub function
function decodingAcc = calc_decoding_acc(idx,K,Events,Time)
T = length(Events);
decodingAcc = zeros(19,1);
for Delay = 0:18
    
    Time_Delay = Time-Delay;
    idx0 = idx(length(find(Time_Delay<1))+1:end-length(find(Time_Delay>T)));% cutoff index of clustering results
    
    Time_Delay(Time_Delay<1)=[];Time_Delay(Time_Delay>T)=[];
    Events_Delay = Events(Time_Delay);% corresponding condition
    
    predict_Events = zeros(size(Events_Delay));
    
    for k = 1:K
        predict_Events(idx0==k) = mode(Events_Delay(idx0==k));
    end
    decodingAcc(Delay+1)=sum(predict_Events==Events_Delay)./length(Events_Delay);
end
end