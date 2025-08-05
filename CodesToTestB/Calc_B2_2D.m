% Calc_B2.m calculates the Jacobian matrix B2 (Algorithm 3 in Vatankhah et al., 2022).
%
%
% Based on Saeed Vatankhah, Institute of Geophysics and Geomatics, China University of Geosciences (Wuhan), October, 2019.
% Udpated Renaut 2022, 2025 for 2D ( But looks like an original error)
%% Input Parameters
% dx: Dimension of prisms in East direction (in meter).
% dy: Dimension of prisms in North direction (in meter).
% nx: Total number of prisms in East direction.
% ny: Total number of prisms in North direction.
% n: Total number of prisms.
% Calculate for each dimension - here for susceptibility
function [B2]=Calc_B2_2D(dx,dy,nx,ny, dxm1,dym1)
dxm1=reshape(dxm1,nx,ny);
dym1=reshape(dym1,nx,ny);
n=nx*ny;
B2=sparse(zeros(n));
s=1;
for j=1:ny %Diagonal block
    i_index=s:s+nx-1;
    B2(i_index,i_index)=sparse(diag(dym1(:,j)/dx-dxm1(:,j)/dy));
    % Next block j+1
    if j<ny
        j_index=s+nx:s+2*nx-1;
        B2(i_index,j_index)=sparse(diag(dxm1(:,j)/dy));
    end
    % i+1 block
    B2(i_index,i_index)=B2(i_index,i_index)-sparse(diag(dym1(1:nx-1,j)/dx,1));
    s=s+nx;
end
end

