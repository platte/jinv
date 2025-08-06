%% Test cross-product linearization

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
[dxm1,dym1,dxm2,dym2] = Deval(D,x,n);

% Calculate cross-product:
tx=Calc_tvec_2D(dxm1,dym1,dxm2,dym2);

% Cross product Jacobian
B1=Calc_B1_2D(hx,hy, nx,ny, dxm2,dym2);
B2=Calc_B2_2D(hx,hy, nx,ny, dxm1,dym1);

%% Generate perturbed x
EP = 10.^(-9:.5:0); % relative perturbation
err = zeros(size(EP)); % error in linearization

for k=1:length(EP)
    
    % generate perturbation in x (call it xp)
    ep = EP(k);
    xr = rand(size(x));
    dx= ep*xr/norm(xr);
    xp = x+dx; % perturbed x

    % Calculate the derivatives for the perturbed state xp
    [dxp1,dyp1,dxp2,dyp2] = Deval(D,xp,n);

    % Calculate the derivatives for the perturbation dx
    [p1x,p1y,p2x,p2y] = Deval(D,dx,n);

    % Calculate cross-product:
    txp=Calc_tvec_2D(dxp1,dyp1,dxp2,dyp2);

    % Linear Approximation of the cross-product using Jacobian
    Bdx = dym2.*p1x-dxm2.*p1y+dxm1.*p2y-dym1.*p2x;
    txpLinear = tx+Bdx; 
   
    err(k) = norm(txpLinear - txp)/norm(txp);
end
figure(100)
loglog(EP,err,'r*-',EP,EP.^2,'--')
xlabel('Perturbation','Interpreter','latex');
ylabel('Error in Linearization','Interpreter','latex');
title('Error vs Perturbation','Interpreter','latex');
legend('error', '$O(\Delta x^2)$','Interpreter','latex')
set(gca,'FontSize',14)
grid on;