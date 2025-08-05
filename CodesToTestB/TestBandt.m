% a basic code to set up the matrices for the matrices B and the cg t
% Rosemary Renaut June 8 2025
% Replace 3D by 2D versions
% hx, hy, grid in x and y directions
% hz=1 for 2D
% nx, ny number of points in grid x and y directions, nz=1 - no z dimension
% n=nx*ny - total number of unknowns
% dim=2 - number of dimensions
% g  (xmodel1) - first model
% m (xmodel2)  - second model
% Dg, Dm - derivatives of model 1 and model 2
% scaling if needed for model 1 and model 2 
scaling=[1,1];
nz=1;dim=2;
nx=5;ny=5;hx=1/(nx-1);hy=1/(ny-1);
n=nx*ny;
Truegx=rand(nx,nx);
Truemx=rand(nx,nx);
x=[Truegx;Truemx];
 Dg = dsOperator('finite difference',  [nx,ny,nz],1);% 
 Dm = dsOperator('finite difference', [nx,ny,nz], 1);
blockmatrix{1}=Dg;blockmatrix{2}=Dm; 
D=BlockMatrixOperator(2,blockmatrix,scaling);% the block D matrix
% Calculate the derivatives in x and y directions for 2 models
Derivx=D*x(:);%assumes that D is square in each dimension
% My original code has square derivatives
ind_d0=1; ind_d1=n;
dxm1=Derivx( ind_d0:ind_d1,:);
ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
dym1=Derivx(  ind_d0:ind_d1,:);
ind_d1=dim*n;
ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
dxm2=Derivx( ind_d0:ind_d1,:);
ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
dym2=Derivx( ind_d0:ind_d1,:);
B1=Calc_B1_2D(hx,hy, nx,ny, dxm2,dym2);
B2=Calc_B2_2D(hx,hy, nx,ny, dxm1,dym1);
figure, tiledlayout('flow')
for k=1:nx
j=(k-1)*nx+1;
nexttile,
spy(B1(j:j+nx-1,j:j+nx-1))
end
figure, tiledlayout('flow')
for k=1:nx
j=(k-1)*nx+1;
nexttile,
spy(B2(j:j+nx-1,j:j+nx-1))
end
B=[B1;B2];
t=Calc_tvec_2D(dxm1,dym1,dxm2,dym2);
tnorm=norm(t)^2;
