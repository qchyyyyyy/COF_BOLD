function [mask,label,oi_x,oi_y,V_mask] = load_cortical_atlas(parcellation,hemisphere,dataRoot)
%MASKLABEL Build a cortical mask and parcel-label outline.
%
% The optional dataRoot argument makes template lookup independent of the
% current MATLAB working directory. Legacy two-argument calls continue to
% use the repository-relative data/ directory.

if nargin < 3 || isempty(dataRoot)
    dataRoot = getenv('COF_DATA_ROOT');
    if isempty(dataRoot), dataRoot='data'; end
end
if isstring(dataRoot)
    dataRoot = char(dataRoot);
end
if ~ischar(dataRoot) || ~isrow(dataRoot)
    error('dataRoot must be a character vector or scalar string.');
end

templateRoot=getenv('COF_TEMPLATE_ROOT');
if isempty(templateRoot), templateRoot=fullfile(dataRoot,'templates'); end
labelIdFile = fullfile(templateRoot, 'HCPex_v1.0', 'HCPex_LabelID.mat');
coordinateFile = fullfile(templateRoot, 'coordinate_cii_v2.mat');
yeo7File = fullfile(templateRoot, 'Yeo2011', 'Yeo2011_7networks.mat');
yeo17File = fullfile(templateRoot, 'Yeo2011', 'Yeo2011_17networks.mat');

if strcmp(hemisphere,'lh')
    hemisphere='L';
end
if strcmp(hemisphere,'rh')
    hemisphere='R';
end
assert(ismember(string(parcellation),["HCPex_360","HCPex_22","Yeo2011_7networks","Yeo2011_17networks"]),'cof:Atlas','Unknown atlas.');
assert(ismember(string(hemisphere),["L","R"]),'cof:Hemisphere','Unknown hemisphere.');
load(labelIdFile,'LabelID')
switch parcellation
    case 'HCPex_360'
        switch hemisphere
            case 'L'
                load(coordinateFile,'cii_mask_left','cii_label_left');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load(coordinateFile,'cii_mask_right','cii_label_right');
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
    case 'HCPex_22'
        switch hemisphere
            case 'L'
                load(coordinateFile,'cii_mask_left','cii_label_left');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load(coordinateFile,'cii_mask_right','cii_label_right');
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
        [~,index] = sort(cell2mat(LabelID(:,2)));
        LabelID_sorted = LabelID(index,:);
        label_22 = zeros(size(label));[M,N] = size(label_22);
        for i = 1:M
            for j = 1:N
                if label(i,j)~=0
                    label_22(i,j) = LabelID_sorted{label(i,j),5};
                end
            end
        end
        label = label_22;
    case 'Yeo2011_7networks'
        switch hemisphere
            case 'L'
                load(coordinateFile,'cii_mask_left','cii_label_left');
                load(yeo7File);
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load(coordinateFile,'cii_mask_right','cii_label_right');
                load(yeo7File);
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
        case 'Yeo2011_17networks'
        switch hemisphere
            case 'L'
                load(coordinateFile,'cii_mask_left','cii_label_left');
                load(yeo17File);
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load(coordinateFile,'cii_mask_right','cii_label_right');
                load(yeo17File);
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
end
outline = conv2(label,ones(2)/4,'valid');
outline_ind = find(round(outline,4)~=label(1:end-1,1:end-1));
[oi_y,oi_x] = ind2sub(size(outline),outline_ind);
oi_y = oi_y+0.5;oi_x = oi_x+0.5;

V_mask = conv2(mask,ones(2,2)/4,'valid');
V_mask(V_mask~=1)=0;

end
