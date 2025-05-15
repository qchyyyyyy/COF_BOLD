clear;
addpath(genpath('functions'));
DataPath = 'data/fMRI/rfMRI_REST1_LR/';
OutputPath = 'data/results/decomposition_results/';

parcellation = 'Yeo2011_7networks';%'Yeo2011_7networks''HCPex_22''Yeo2011_17networks''HCPex_360'
hemisphere = 'LR';
[mask_L,label_L,oi_x_L,oi_y_L,V_mask_L] = masklabel(TemplatePath,parcellation,hemisphere(1));
[mask_R,label_R,oi_x_R,oi_y_R,V_mask_R] = masklabel(TemplatePath,parcellation,hemisphere(2));

mask = [mask_L mask_R];
label = [label_L label_R];
oi_x = [oi_x_L;oi_x_R+251];
oi_y = [oi_y_L;oi_y_R];
V_mask = [V_mask_L V_mask_R];
clearvars mask_L mask_R label_L label_R oi_x_L oi_x_R oi_y_L oi_y_R V_mask_L V_mask_R
%%

File = dir(fullfile(DataPath));
SubjectList = cell(size(File));
for m = 1:length(File)
    SubjectList{m}=File(m).name;
end
SubjectList(1:2)=[];

%%
Nodes = nnz(find(V_mask==1));
AAT = zeros(Nodes,Nodes);

for sub = length(SubjectList):-1:1
    disp(['sub = ' num2str(sub)]);
    patternProfile = generate_patternProfile(...
        [DataPath,SubjectList{sub},'/'],...
        V_mask); % generate Nodes*T matrix
    AAT = AAT + patternProfile*patternProfile';
end


save([OutputPath hemisphere '_normDivSIGN_AAT_group_REST.mat'],'AAT','-v7.3');
%%
[U,S2] = eig(AAT);
S = sqrt(diag(S2));clearvars S2

save([OutputPath hemisphere '_normDivSIGNMode_group_REST_LR.mat'],...
    'U','S','-v7.3');
%%

function DivCurl = generate_patternProfile(Path,V_mask)

L = load([Path,'VelField_lh_hilbert.mat'],'VelField');
R = load([Path,'VelField_rh_hilbert.mat'],'VelField');

Start = 101;
T = min(1000,size(L.VelField.Ux,3)-Start-98);

load([Path,'TC_mesh_r_f_zsore_s.mat'],'data');
BrainImg = cat(2,data.BrainImg_lh(:,1:end-1,Start:Start+T),data.BrainImg_rh(:,:,Start:Start+T));clearvars data
BrainImg = convn(BrainImg,ones(2,2,2)/8,'valid');
SignBOLD = sign(BrainImg);
SignBOLD = reshape(SignBOLD,[],T);
SignBOLD = SignBOLD(V_mask==1,:);

filtSigma = 3;
filtWidth = round(3*filtSigma);
imageFilter=fspecial('gaussian',filtWidth,filtSigma);

DivMap = [];
Ux = cat(2,L.VelField.Ux(:,:,Start:Start+T-1),R.VelField.Ux(:,:,Start:Start+T-1));
Uy = cat(2,L.VelField.Uy(:,:,Start:Start+T-1),R.VelField.Uy(:,:,Start:Start+T-1));clearvars L R
for t = T:-1:1
    div_oneframe= divergence(Ux(:,:,t),Uy(:,:,t));
    div_oneframe(V_mask==0)=nan;
    DivMap(:,:,t) = nanconv(div_oneframe,imageFilter);
end
clearvars Ux Uy
DivMap(repmat(V_mask,1,1,T)==0)=0;

DivMap = reshape(DivMap,[],T);
DivMap = DivMap(V_mask==1,:);

DivCurl = zscore(DivMap,[],2) + 1i*SignBOLD;
end