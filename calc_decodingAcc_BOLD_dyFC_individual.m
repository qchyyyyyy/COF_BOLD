clear;
addpath(genpath('functions'));

Condition = {'WM','RELATIONAL'};
PE = {'LR'};
parcellation = 'HCPex_360';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';

[mask_L, label_L, oi_x_L, oi_y_L, V_mask_L] = masklabel(parcellation,hemisphere(1));
[mask_R, label_R, oi_x_R, oi_y_R, V_mask_R] = masklabel(parcellation,hemisphere(2));

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R

label(isnan(label))=0;
label_vec = label(mask==1);


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
                T = 405;
            case 'RELATIONAL'
                T = 232;
        end
        SubjectList(BadInd)=[];
        
        Events = zeros(T,1);
        Time = 10:T-9;
        
        for bb = 1:8
            idx = ceil(meanStimType(bb,1)):fix(meanStimType(bb,2));
            idx(idx>T)=[];idx(idx<1)=[];
            Events(idx) = meanStimType(bb,3);
        end
        
        Events(Events>0 & Events<5)=1;
        Events(Events>4)=2;
        Events = Events+1;
        
        
        BrainImg_group=zeros([size(mask),T]);
%         dyFC_group = 
        
        for sub = length(SubjectList):-1:1
            disp(['sub = ' num2str(sub)]);
            load([DataPath,SubjectList{sub}, '/TC_mesh_r_f_zscore.mat']);
            BrainImg = cat(2,data.BrainImg_lh,data.BrainImg_rh);clearvars data
            

            Feature = map2vec(BrainImg(:,:,Time),mask);
            Feature = zscore(Feature,[],2);Cov = corr(Feature,'rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_BOLD(:,sub) = calc_decoding_acc(idx,K,Events,Time);
            
            BrainImg = map2vec(BrainImg,mask);
            for i_re = 1:max(label,[],'all')
                 BrainImg_region(i_re,:) = mean(BrainImg(label_vec==i_re,:));
            end
            Feature = coupling(BrainImg_region',14);
            FC_mask = ones(max(label,[],'all'),max(label,[],'all'));
            FC_mask = triu(FC_mask,1);
            Feature = map2vec(Feature,FC_mask);
            
            Cov = corr(Feature(:,Time),'rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_dyFC(:,sub) = calc_decoding_acc(idx,K,Events,Time);
            
            Feature = zscore(Feature,[],2);
            Cov = corr(Feature(:,Time),'rows','pairwise');
            rng(0);idx = community_louvain(Cov,1,[],'negative_asym');K = max(idx);
            decodingAcc_dyFC_zscore(:,sub) = calc_decoding_acc(idx,K,Events,Time);
        end


        save([OutputPath 'decoding_results/' hemisphere '_decodingAcc_BOLD_individual_' Condition{cdn} '_' PE{pe} '.mat'],...
            'decodingAcc_BOLD','-v7.3');
        save([OutputPath 'decoding_results/' hemisphere '_decodingAcc_dyFC_individual_' Condition{cdn} '_' PE{pe} '.mat'],...
            'decodingAcc_dyFC','decodingAcc_dyFC_zscore','-v7.3');
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