clear;
N = 17;
Yeo2011 = ft_read_cifti(['D:\spatiotemporal patterns\cifti data\Yeo2011\Yeo2011_',num2str(N),'Networks_N1000.dlabel.nii']);
label = Yeo2011.parcels;
% ROIname = Yeo2011.parcelslabel;

load('D:\spatiotemporal patterns\spatiotemporal pattrens\DataPreprocessing\coordinate_cii_v2.mat');
mask_region_right = label(Yeo2011.brainstructure==2);mask_region_right(i_invalid_right)=[];
mask_region_left = label(Yeo2011.brainstructure==1);mask_region_left(i_invalid_left)=[];

downSRate = 2;
x_target_cood = -250:downSRate:250;y_target_cood = -200:downSRate:200;

cii_label_right = nan(length(y_target_cood),length(x_target_cood));
for i_x = 1:length(x_target_cood)
    for i_y = 1:length(y_target_cood)
        for i_point = 1:length(x_right)
            d(i_point) = sqrt((x_target_cood(i_x)-x_right(i_point)).^2+(y_target_cood(i_y)-y_right(i_point)).^2);
        end
        [~,i_roi] = min(d);
        cii_label_right(i_y,i_x) = mask_region_right(i_roi);
    end
end
cii_label_right = cii_label_right.*cii_mask_right;

cii_label_left = nan(length(y_target_cood),length(x_target_cood));
for i_x = 1:length(x_target_cood)
    for i_y = 1:length(y_target_cood)
        for i_point = 1:length(x_left)
            d(i_point) = sqrt((x_target_cood(i_x)-x_left(i_point)).^2+(y_target_cood(i_y)-y_left(i_point)).^2);
        end
        [~,i_roi] = min(d);
        cii_label_left(i_y,i_x) = mask_region_left(i_roi);
    end
end
cii_label_left = cii_label_left.*cii_mask_left;
save(['D:\spatiotemporal patterns\cifti data\Yeo2011\Yeo2011_',num2str(N),'networks.mat'],'cii_label_left','cii_label_right');