% Calc_tvec.m produces cross-gradient vector "t" 
%
% Based Originally on Code Written by Saeed Vatankhah, Institute of Geophysics and Geomatics, China University of Geosciences (Wuhan), October, 2019.  
% Updated to just use derivatives Renaut 2022 June 8 2025 
%% Input Parameters
function [t]=Calc_tvec_2D(dxm1,dym1,dxm2,dym2)
% Ignoring scaling here since both terms have the same denominator and
% scaling is absorbed into regularization parameter
t=dxm1(:).*dym2(:)-dym1(:).*dxm2(:);
end

