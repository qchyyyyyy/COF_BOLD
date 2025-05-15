function [mask,label,oi_x,oi_y,V_mask] = masklabel(parcellation,hemisphere)
if strcmp(hemisphere,'lh')
    hemisphere='L';
end
if strcmp(hemisphere,'rh')
    hemisphere='R';
end
load('data/templates/HCPex_v1.0/HCPex_LabelID.mat')
switch parcellation
    case 'HCPex_360'
        switch hemisphere
            case 'L'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_label_left');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_right','cii_label_right');
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
    case 'HCPex_22'
        switch hemisphere
            case 'L'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_label_left');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_right','cii_label_right');
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
                load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_label_left');
                load('data/templates/Yeo2011/Yeo2011_7networks.mat');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_right','cii_label_right');
                load('data/templates/Yeo2011/Yeo2011_7networks.mat');
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
        case 'Yeo2011_17networks'
        switch hemisphere
            case 'L'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_left','cii_label_left');
                load('data/templates/Yeo2011/Yeo2011_17networks.mat');
                mask = cii_mask_left;label = cii_label_left;
            case 'R'
                load('data/templates/coordinate_cii_v2.mat','cii_mask_right','cii_label_right');
                load('data/templates/Yeo2011/Yeo2011_17networks.mat');
                mask = cii_mask_right;label = cii_label_right;
        end
        label(isnan(label))=0;
        mask(isnan(mask))=0;
end
outline = conv2(label,ones(2)/4,'valid');
outline_ind = find(round(outline,4)~=label(1:end-1,1:end-1));
[oi_y,oi_x] = ind2sub(size(outline),outline_ind);
oi_y = oi_y+0.5;oi_x = oi_x+0.5;

V_mask = zeros(size(mask)-1);
[M,N] = size(V_mask);
for i = 1:M
    for j = 1:N
        if mask(i,j)==1 && mask(i,j+1)==1 && mask(i+1,j)==1 && mask(i+1,j+1)==1
            V_mask(i,j) = 1;
        end
    end
end
end