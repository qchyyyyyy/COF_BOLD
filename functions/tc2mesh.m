function TC_mesh = tc2mesh(TC,x,y,mask)
x = double(x);y = double(y);
downSRate = 2 ;
[Ylength,Xlength] = size(mask);
if Ylength==502
    x_target_cood = 2-Xlength:downSRate:Xlength-2;
else
    x_target_cood = 1-Xlength:downSRate:Xlength-1;
end
y_target_cood = 1-Ylength:downSRate:Ylength-1;
[xi,yi] = meshgrid(x_target_cood,y_target_cood);

for i_timepoint = size(TC,2):-1:1
    F_TC = scatteredInterpolant(x,y,TC(:,i_timepoint));
    TC_mesh(:,:,i_timepoint) = F_TC(xi, yi).*mask;
end