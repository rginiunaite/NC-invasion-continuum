function TwoD_ComputeWeightsUsingTwoDifferentMethods

%{
 diary('TwoD_ComputeWeightsUsingTwoDifferentMethods.diary');
 TwoD_ComputeWeightsUsingTwoDifferentMethods
 diary off
%}

oldpath = path;
addpath('..');

R=1.1;
h=1/25;
N1=70;
N2=90;
NR=4;
tol=1e-2;
OmegaStr = '@(r)(3/(pi*R^2)*(1-r/R))'; 

disp('WARNING: This function may take a lot of time to run!!!');
disp(' ');

disp('Using composite trapezoidal rule for x1-direction ...')
x1DirTrap = setupIntegralRule2D('x1', R, h, N1, N2, eval(OmegaStr), ...
				      1, NR, tol,true);
disp('Using composite trapezoidal rule for x2-direction ...')
x2DirTrap = setupIntegralRule2D('x2', R, h, N1, N2, eval(OmegaStr), ...
				      1, x1DirTrap.NR, ...
				      inf, true);


disp('Using composite midpoint rule for x2-direction ...')
x2DirMid = setupIntegralRule2D('x2', R, h, N1, N2, eval(OmegaStr), ...
				     2, NR, tol,true);
disp('Using composite midpoint rule for x1-direction ...')
x1DirMid = setupIntegralRule2D('x1', R, h, N1, N2, eval(OmegaStr), ...
				     2, x2DirMid.NR, ...
				     inf, true);
rmpath('..')