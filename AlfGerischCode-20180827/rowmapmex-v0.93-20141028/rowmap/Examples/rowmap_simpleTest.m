function rowmap_simpleTest
% function rowmap_simpleTest
%
% This example script solves a simple (diagonal) 3-dimensional
% constant coefficient linear ODE system with the ROWMAP time
% integration scheme for various tolerances. An example of an
% output function is given, which serves as a basis for
% user-supplied, more advanced output functions. These might
% store and/or post-process particular solution values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : rowmap_simpleTest.m
% Version : 25 January 2009 (Alf Gerisch, University of Halle)
%
% Documentation of the ROWMAP MEX Interface is maintained at
% http://sim.mathematik.uni-halle.de/~gerisch/r/rowmapmex.html
%
% Please send comments and bug reports to 
%     alf.gerisch@mathematik.uni-halle.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format compact
format short e

addpath('..'); % path to rowmap

%
% set tspan vector for call of time integrator
%
tspan = [0:0.2:1];

%
% set initial value for ODE system
%
y0    = ones(3,1);

%
% set options for rowmap call
%
clear options;
options.OptWarnMiss = 0;       % set to 1 to see more possible options
options.RHSnonautonomous = 0;  % the ODE system is autonomous

% set some options for the time step selection etc. (others like
% tolerances are in the loop below)
options.MaxNumberOfSteps=1000;
options.InitialStep=1e-2;

% select what we want to have returned in t and y
options.ReturnMode = 0;  % 0 .. empty t and y returned
options.ReturnMode = 1;  % 1 .. t_end and corresponding approximation

% set output function and when it is to be called
options.OutputFun = @localOutputfun; % is defined below
options.OutputCallMode = 4; % 4 .. call after each time step and for
                            % each tspan value
%options.OutputCallMode = 0; % 0 .. call never

% set post-step function and when it is to be called
%options.PostStepFun = @; % we provide no post-step function here 
options.PostStepCallMode = 0; % 0 .. means never

tols = 10.^[-3:-1:-8];
allResults = nan(length(tols),5);
for i=1:length(tols)
  tol = tols(i);
  disp(' ');
  disp('--->');
  disp(['---> starting for tolerance tol = ' num2str(tol)]);
  disp('--->');
  options.RelTol=tol;
  options.AbsTol=tol;

  [t, y, stat, hnext] = rowmap(@rhs, tspan, y0, options);
  %tol, t, y, stat, hnext
  
  sol = y0.*exp([-1:-1:(-length(y0))]'*tspan(end));
  err=norm(y(end,:)'- sol);
  allResults(i,:)=[tol err stat(2) stat(3) stat(1)];
end

disp(' ');
disp('--->');
disp('---> Result table');
disp('--->');
disp('   tol          err          nsteps        nstepsrej   IDID')
disp('-------------------------------------------------------------------')
format short e
disp(allResults)


rmpath('..')

return

% End of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ydot = rhs(t,y)

ydot=[-1:-1:(-length(y))]'.*y;

return
% End of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = localOutputfun(t, y, flag)
% sample output function
%
% this function can be used to store solution values, generate
% animations etc.

status = 0; % means continue time integration
%status = 1; % means return from time integration
  
if (strcmp(flag,'init'))
  disp('localOutputfun::Received flag ''init'' in output function.');
  disp(['   initial time   = ' num2str(t(1))]);
  disp(['   initial vector = ' num2str(y')]);
  return
end

if (strcmp(flag,'done'))
  disp('localOutputfun::Received flag ''done'' in output function.');
  return
end


% otherwise (flag == ''): y is vector with approximate 
% solution at time t  
disp(['localOutputfun::Reached time t = ' num2str(t)]);

return
% End of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
