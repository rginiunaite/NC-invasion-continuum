%function test_fft2

h = 0.1;
R = 1;
N1 = 200;
N2 = 220;
Omega = 1;
ruleId = 2;
NR = 20;

% use old methods to setup
BCs = 'pp';
clear mask_old
mask_old.x1Dir = setupIntegralRule2D_old('x1', R, h, N1, N2, Omega, ...
					       ruleId, NR);
mask_old.x2Dir = setupIntegralRule2D_old('x2', R, h, N1, N2, Omega, ...
					       ruleId, NR);
% use new methods to setup 
BCs = 'pp';
clear mask_pp
mask_pp.BCs = BCs;
nonlocal = setupIntegralRule2D_weights('x1', R, h, Omega, ...
					     ruleId, NR);
mask_pp.x1Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);
nonlocal = setupIntegralRule2D_weights('x2', R, h, Omega, ...
					     ruleId, NR);
mask_pp.x2Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);

% use new methods to setup pppp
BCs = 'pppp';
clear mask_pppp
mask_pppp.BCs = BCs;
nonlocal = setupIntegralRule2D_weights('x1', R, h, Omega, ...
					     ruleId, NR);
mask_pppp.x1Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);
nonlocal = setupIntegralRule2D_weights('x2', R, h, Omega, ...
					     ruleId, NR);
mask_pppp.x2Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);

% use new methods to setup zzzz
BCs = 'zzzz';
clear mask_zzzz
mask_zzzz.BCs = BCs;
nonlocal = setupIntegralRule2D_weights('x1', R, h, Omega, ...
					     ruleId, NR);
mask_zzzz.x1Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);
nonlocal = setupIntegralRule2D_weights('x2', R, h, Omega, ...
					     ruleId, NR);
mask_zzzz.x2Dir = setupIntegralRule2D_BCs(nonlocal, N1, N2, BCs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate integral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = rand(N2, N1);

methodId = 1;
[A11_old_1, A22_old_1] = evalIntegral2D_old(G, mask_old, methodId);
methodId = 2;
[A11_old_2, A22_old_2] = evalIntegral2D_old(G, mask_old, methodId);
methodId = 3;
[A11_old_3, A22_old_3] = evalIntegral2D_old(G, mask_old, methodId);

methodId = 1;
[A11_pp_1, A22_pp_1] = evalIntegral2D(G, mask_pp, methodId);
methodId = 2;
[A11_pp_2, A22_pp_2] = evalIntegral2D(G, mask_pp, methodId);
methodId = 3;
[A11_pp_3, A22_pp_3] = evalIntegral2D(G, mask_pp, methodId);

methodId = 1;
[A11_pppp_1, A22_pppp_1] = evalIntegral2D(G, mask_pppp, methodId);
disp(['max norm periodicity method 1: ' num2str(norm(A11_pppp_1(:,1)-A11_pppp_1(:,end),inf))]);
A11_pppp_1 =A11_pppp_1(:,2:end); 
disp(['max norm periodicity method 1: ' num2str(norm(A22_pppp_1(1,:)-A22_pppp_1(end,:),inf))]);
A22_pppp_1 =A22_pppp_1(2:end,:); 
methodId = 3;
[A11_pppp_3, A22_pppp_3] = evalIntegral2D(G, mask_pppp, methodId);
disp(['max norm periodicity method 3: ' num2str(norm(A11_pppp_3(:,1)-A11_pppp_3(:,end),inf))]);
A11_pppp_3 =A11_pppp_3(:,2:end); 
disp(['max norm periodicity method 3: ' num2str(norm(A22_pppp_3(1,:)-A22_pppp_3(end,:),inf))]);
A22_pppp_3 =A22_pppp_3(2:end,:); 

methodId = 1;
[A11_zzzz_1, A22_zzzz_1] = evalIntegral2D(G, mask_zzzz, methodId);
methodId = 3;
[A11_zzzz_3, A22_zzzz_3] = evalIntegral2D(G, mask_zzzz, methodId);

format compact

disp('A11')
disp('--------------------------------');
ref = A11_old_2;

disp('old pp')
disp([isreal(A11_old_1) norm(imag(A11_old_1(:)),inf) ...
 norm(A11_old_1(:)-ref(:),inf)])
disp([isreal(A11_old_2) norm(imag(A11_old_2(:)),inf) ...
 norm(A11_old_2(:)-ref(:),inf)])
disp([isreal(A11_old_3) norm(imag(A11_old_3(:)),inf) ...
 norm(A11_old_3(:)-ref(:),inf)])

disp('new pp')
disp([isreal(A11_pp_1) norm(imag(A11_pp_1(:)),inf) ...
 norm(A11_pp_1(:)-ref(:),inf)])
disp([isreal(A11_pp_2) norm(imag(A11_pp_2(:)),inf) ...
 norm(A11_pp_2(:)-ref(:),inf)])
disp([isreal(A11_pp_3) norm(imag(A11_pp_3(:)),inf) ...
 norm(A11_pp_3(:)-ref(:),inf)])

disp('new pppp')
disp([isreal(A11_pppp_1) norm(imag(A11_pppp_1(:)),inf) ...
 norm(A11_pppp_1(:)-ref(:),inf)])
disp([isreal(A11_pppp_3) norm(imag(A11_pppp_3(:)),inf) ...
 norm(A11_pppp_3(:)-ref(:),inf)])


disp('')

disp('A22')
disp('--------------------------------');
ref = A22_old_2;

disp('old pp')
disp([isreal(A22_old_1) norm(imag(A22_old_1(:)),inf) ...
 norm(A22_old_1(:)-ref(:),inf)])
disp([isreal(A22_old_2) norm(imag(A22_old_2(:)),inf) ...
 norm(A22_old_2(:)-ref(:),inf)])
disp([isreal(A22_old_3) norm(imag(A22_old_3(:)),inf) ...
 norm(A22_old_3(:)-ref(:),inf)])

disp('new pp')
disp([isreal(A22_pp_1) norm(imag(A22_pp_1(:)),inf) ...
 norm(A22_pp_1(:)-ref(:),inf)])
disp([isreal(A22_pp_2) norm(imag(A22_pp_2(:)),inf) ...
 norm(A22_pp_2(:)-ref(:),inf)])
disp([isreal(A22_pp_3) norm(imag(A22_pp_3(:)),inf) ...
 norm(A22_pp_3(:)-ref(:),inf)])

disp('new pppp')
disp([isreal(A22_pppp_1) norm(imag(A22_pppp_1(:)),inf) ...
 norm(A22_pppp_1(:)-ref(:),inf)])
disp([isreal(A22_pppp_3) norm(imag(A22_pppp_3(:)),inf) ...
 norm(A22_pppp_3(:)-ref(:),inf)])

disp('')

disp('A11 zzzz')
disp('--------------------------------');
ref = A11_zzzz_1;
disp([isreal(A11_zzzz_1) norm(imag(A11_zzzz_1(:)),inf) ...
 norm(A11_zzzz_1(:)-ref(:),inf)])
disp([isreal(A11_zzzz_3) norm(imag(A11_zzzz_3(:)),inf) ...
 norm(A11_zzzz_3(:)-ref(:),inf)])

disp('A22 zzzz')
disp('--------------------------------');
ref = A22_zzzz_1;
disp([isreal(A22_zzzz_1) norm(imag(A22_zzzz_1(:)),inf) ...
 norm(A22_zzzz_1(:)-ref(:),inf)])
disp([isreal(A22_zzzz_3) norm(imag(A22_zzzz_3(:)),inf) ...
 norm(A22_zzzz_3(:)-ref(:),inf)])


