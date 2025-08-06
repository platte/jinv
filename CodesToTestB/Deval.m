function [dxm1,dym1,dxm2,dym2] = Deval(D,x,n)

Derivx=D*x(:);%assumes that D is square in each dimension
% My original code has square derivatives
ind_d0=1; ind_d1=n;
dxm1=Derivx( ind_d0:ind_d1,:);
ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
dym1=Derivx(  ind_d0:ind_d1,:);
ind_d1=2*n;
ind_d0=1+ind_d1;  ind_d1=ind_d1+n;
dxm2=Derivx( ind_d0:ind_d1,:);
ind_d0=ind_d1+1;  ind_d1=ind_d1+n;
dym2=Derivx( ind_d0:ind_d1,:);

end