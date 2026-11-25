function result = wrap_angle_difference(a,b,angular)
% Signed elementwise difference, optionally mapped to [-pi, pi).
% Release implementation of elementary periodic arithmetic (MIT).
if nargin<3, angular=true; end
result=a-b;
if angular==1, result=mod(result+pi,2*pi)-pi; end
end
