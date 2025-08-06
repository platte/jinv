function checkDerivative(fcn, x, secondDerivative)
% checkDerivative - numerical test of derivative
%
% function checkDerivative(fcn, x);
%
% Input arguments:
%   fcn   - function
%   x     - evaluation point
%
% Description
%   Test derivative via comparison of linear and quadratic Taylor approximation

if nargin < 3
  secondDerivative = 0;
end

% close all
fprintf('Check Derivative via comparison of Taylor series T0 and T1\n\n');
fprintf('T0 = |f - fh|    and    T1 = |f + h * df  - fh|\n\n');
fprintf('    h              T0            T1 \n');

% initial function and gradient/Jacobian evaluation at x
if secondDerivative
  y = randn(size(x));
  fcn = @(x) checkSecondDerivative(fcn,x,y);
end

[f, df] = feval(fcn, x); f = f(:);

% evaluation grid
h = logspace(1,-10,20);

% choose random direction
v = randn(size(x)); T0 = zeros(size(h)); T1 = T0;
 
for j=1:length(h)
	
	% function evaluation at x+h(j)*v
  fh = feval(fcn, x + h(j)*v ); fh = fh(:);
	% Taylor approximation of first and second order
  T0(j) = norm(f - fh); T1(j) = norm(f + h(j) * df' * v - fh);

  fprintf('%10.4e     %10.4e    %10.4e\n', h(j), T0(j), T1(j));
end;

% plot results
figure,
loglog(h, T0, 'b-'), hold on
loglog(h, T1, 'r--'),
title(['Taylor approximation $|T_0(h) - f|$ and $|T_1(h)  - f|$'],'Interpreter','Latex')
legend('T_0','T_1','Location','northwest')
xlabel('$h$'), ylabel('$T_0(h) = |f-fh|$, $T_1(h) = |f + h f'' - fh|$','Interpreter','Latex')
% text(func2str(fcn))
text(0.25,0.95,['f = ', func2str(fcn)],'Units','normalized')
drawnow

return
 
function [f, df] = checkSecondDerivative(fcn,x,y)
  [~, g, H] = fcn(x);
  f = g'*y; df = H*y;
return
