% Calc_B1.m calculates the Jacobian matrix B1 (Algorithm 3 in Vatankhah et al., 2022).
%
% Based on Saeed Vatankhah, Institute of Geophysics and Geomatics, China University of Geosciences (Wuhan), October, 2019.
% Udpated Renaut 2022, 2025 for 2D ( But looks like an original error)
%% Input Parameters
% dx: Dimension of prisms in East direction (in meter).
% dy: Dimension of prisms in North direction (in meter).
% nx: Total number of prisms in East direction.
% ny: Total number of prisms in North direction.
% n: Total number of prisms.
% Calculate for each dimension - here for Density
function [B1]=Calc_B1_2D(dx,dy,nx,ny,dxm2,dym2)
dxm2=reshape(dxm2,nx,ny);
dym2=reshape(dym2,nx,ny);
n=nx*ny;
B1=sparse(zeros(n));
s=1;
for j=1:ny
    %Diagonal block
    i_index=s:s+nx-1;
    B1(i_index,i_index)=sparse(diag(dxm2(:,j)/dy-dym2(:,j)/dx));
    % Next block j+1
    if j<ny
        j_index=s+nx:s+2*nx-1;
        B1(i_index,j_index)=-sparse(diag(dxm2(:,j)/dy));
    end
    % i+1 block
    B1(i_index,i_index)=B1(i_index,i_index)+sparse(diag(dym2(1:nx-1,j)/dx,1));
    s=s+nx;
end
end

