function TC_mesh = interpolate_to_grid(TC,x,y,x_target_cood,y_target_cood,mask,method)

if ~exist('mask','var') || isempty(mask)
    mask = ones(length(y_target_cood),length(x_target_cood));
end
if ~exist('method','var') || isempty(method)
    method = 'linear';
end

TC = double(TC);
x = double(x);y = double(y);
x_target_cood = double(x_target_cood);y_target_cood = double(y_target_cood);

[xi,yi] = meshgrid(x_target_cood,y_target_cood);
for i_timepoint = size(TC,2):-1:1
    F_TC = scatteredInterpolant(x,y,double(TC(:,i_timepoint)),method);
    TC_mesh(:,:,i_timepoint) = F_TC(xi, yi).*mask;
end