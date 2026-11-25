function [Ux,Uy] = phase_gradient(phaseSig)
% PHASEGRADIENT calculate the velocity vector fields (VVFs) via phase
% gradient

Ux = zeros(size(phaseSig));
Uy = zeros(size(phaseSig));
for iTime = 1:size(phaseSig,3)
    for iX = 1:size(phaseSig,1)
        Ux(iX,2:end-1,iTime) = -(wrap_angle_difference(phaseSig(iX,3:end,iTime),phaseSig(iX,1:end-2,iTime)))/2;
    end
    for iY = 1:size(phaseSig,2)
        Uy(2:end-1,iY,iTime) = -(wrap_angle_difference(phaseSig(3:end,iY,iTime),phaseSig(1:end-2,iY,iTime)))/2;
    end
end