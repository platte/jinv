%% Test cross-product linearization

%% This is the code from Rosie's TestBandt.m

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
ind_d0=ind_d1+1;  ind_d1=ind_d1+n;
dym2=Derivx( ind_d0:ind_d1,:);

% Calculate cross-product:
tx =  dxm1.*dym2 - dym1.*dxm2;

% Cross product Jacobian
B1=Calc_B1_2D(hx,hy, nx,ny, dxm2,dym2);
B2=Calc_B2_2D(hx,hy, nx,ny, dxm1,dym1);

%% Generate perturbed x
EP = 10.^(-8:.5:0); % relative perturbation
err = zeros(size(EP)); % error in linearization

for k=1:length(EP)
    
    % generate perturbation in x (call it xp)
    ep = EP(k);
    xr = rand(size(x));
    dx= ep*xr/norm(xr);
    xp = x+dx; % perturbed x

    % Calculate the derivatives for the perturbed state xp
    % copied from above
    Derivxp = D*xp(:);
    ind_d0 = 1; ind_d1 = n;
    dxp1 = Derivxp(ind_d0:ind_d1,:);
    ind_d0 = 1 + ind_d1; ind_d1 = ind_d1 + n;
    dyp1 = Derivxp(ind_d0:ind_d1,:);
    ind_d1=dim*n;
    ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
    dxp2=Derivxp( ind_d0:ind_d1,:);
    ind_d0=ind_d1+1;  ind_d1=ind_d1+n;
    dyp2=Derivxp( ind_d0:ind_d1,:);

    % Compute the perturbed cross product
    txp = dxp1.*dyp2 - dyp1.*dxp2;

    % Calculate the difference between perturbed and original cross products
    crossProductDifference = txp - tx;
    norm(crossProductDifference);

    txpLinear = tx+[B1 B2]*dx(:);
    err(k) = norm(txpLinear - txp)/norm(txp);
end
figure(100)
loglog(EP,err,'r*-',EP,EP,'--')