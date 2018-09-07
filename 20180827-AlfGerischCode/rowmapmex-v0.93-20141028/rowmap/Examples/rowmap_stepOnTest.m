function rowmap_stepOnTest
% function rowmap_stepOnTest
%
% This example script solves a simple (diagonal) 3-dimensional
% constant coefficient linear ODE system with the ROWMAP time
% integration scheme AND exactly stepping on some predefined output
% time points. A post-step function is used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : rowmap_stepOnTest.m
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
%options.OutputFun = @; % we do not provide an output function
options.OutputCallMode = 0; % 0 .. call never (just to be sure...)

% set post-step function and when it is to be called
options.PostStepFun = @localPostStepFun; % defined below
options.PostStepCallMode = 0; % 0 .. means call never
options.PostStepCallMode = 3; % 3 .. means call after each successful
                              % time step

tols = 10.^[-6];
allResults = nan(length(tols),5);
for i=1:length(tols)
  tol = tols(i);
  disp(' ');
  disp('--->');
  disp(['---> starting for tolerance tol = ' num2str(tol)]);
  disp('--->');
  options.RelTol=tol;
  options.AbsTol=tol;
  
  global sol
  clear global sol % make sure that there is no global variable sol
                   % before integrations starts; sol is created in
                   % the post-step function.  
  [t, y, stat, hnext] = rowmapStepOn(@rhs, tspan, y0, options);
  %tol, t, y, stat, hnext
  [t',y]
  
  % We could also access the global sol struct generated in the
  % post-step function but we don't do that here.
  
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
function [status, ynew, fnew] = localPostStepFun(t, y, f, flag)
% sample output function
%
% this function can be used to store solution values, generate
% animations etc.

%   There must always be three return values (can be empty). In calls
%   with flag == 'init' or == 'done' all return values are ignored. 
%   In calls with flag == '' the interpretation of ynew and fnew 
%   depends on the value of status.
%   status == 1  ..  abort execution of rowmap
%   status ~= 1  ..  continue with time integration AND
%   if
%   status == 0  ..  no changes to y: ignore ynew and fnew
%   status == 2  ..  changes to y: use ynew to continue integration;
%   	       	      ignore fnew (the new function value is computed 
%	              by an additional call to the user-supplied
% 		      right-hand side function.
%   status == 3  ..  changes to y: use ynew and fnew to continue 
%   	       	      integration; no additional call to the
% 		      user-supplied right-hand side function.

% set default return values
status = 0; 
ynew = [];
fnew = [];

% Store output data in global variable sol. When the computation
% is completed, then the solution can be accessed via the global
% sol variable. 
global sol

if (strcmp(flag,'init'))
  disp('localPostStepFun::Received flag ''init'' in output function.');
  disp(['   initial time   = ' num2str(t(1))]);
  disp(['   initial vector = ' num2str(y')]);
  
  % Store output data in global variable sol. When the computation
  % is completed, then the solution can be accessed via the global
  % sol variable. 
  if ~isfield(sol, 'out')  
    % if sol.out does not exist yet define output time points
    sol.out.tspan = [0.2:0.2:0.7 0.75]; 
    sol.out.count = 0;
    sol.out.t =[];
    sol.out.y =[];
  end
  next = sol.out.count+1;
  if ( next <= length(sol.out.tspan))
    if (abs(sol.out.tspan(next) - t(1)) <= eps ) 
      % save initial data
        sol.out.t(next)   = t(1);
	sol.out.y(next,:) = y';
	sol.out.count = next;
    end
  end
  
  sol.told  = t;   % last time point for interpolation 
  sol.yold  = y;   % last solution for interpolation
  sol.fold  = rhs(t,y); % last rhs side value for cubic interpolation
  return
end

if (strcmp(flag,'done'))
  disp('localPostStepFun::Received flag ''done'' in output function.');

  %disp('...output of [sol.t sol.ypre sol.ypost]');
  %disp([sol.t'  sol.ypre sol.ypost]);
  disp('...output of [sol.out.t sol.out.y]');
  disp([sol.out.t' sol.out.y]);
  
  return
end


% otherwise (flag == ''): y is vector with approximate 
% solution at time t  
disp(['localPostStepFun::Reached time t = ' num2str(t)]);

% here you can modify the y value with which the integration
% continues (set status to 2 or 3 to indicate that) but we do not
% change anything, so
status = 0;

% compute output solution
while (true)
  next = sol.out.count + 1;
  if (next > length(sol.out.tspan))
    break % no more output desired
  end
  if (abs(sol.out.tspan(next) - t) <= eps)
    % we have reached an output time point exactly
    sol.out.t(next)   = t;
    disp(['...storing COMPUTED solution at output time point t = ' ...
	  num2str(sol.out.t(next))]);
    if length(ynew) % here the user can decide what to store --
                    % this here is one possibility
      sol.out.y(next,:) = ynew'; % use ynew already!!
    else
      sol.out.y(next,:) = y'; % use y
    end
    sol.out.count     = next;
    break % all other output times must be in the future
  end
  if (sol.out.tspan(next) > t)
    break % next output time point is still in the future
  end
  % now we need to store solution at told < sol.out.tspan(next) < t
  sol.out.t(next)   = sol.out.tspan(next);
  disp(['...storing INTERPOLATED solution at output time point t = ' ...
	num2str(sol.out.t(next))]);
  % use linear or cubic interpolation (with y and not ynew!!)
  sol.out.y(next,:) = ...
      rowmap_linInterp(sol.told, sol.yold, t, y, sol.out.t(next));
  sol.out.y(next,:) = ...
      rowmap_cubInterp(sol.told, sol.yold, sol.fold, ...
		       t, y, f, sol.out.t(next));
  sol.out.count     = next;
end
% store current time and (modified) solution and rhs value
% as last computed solution
sol.told  = t;      % last time point for interpolation 
if length(ynew) % here the user can decide what to store --
		% this here is the most sensible I think: whenever
                % there is a ynew, we use that in the next interpolation!
  sol.yold = ynew; % use ynew already!!
  sol.fold = fnew; % use fnew already, so it must be computed in
                   % this function (return status = 3) but is only
                   % needed for cibic interpolation
else
  sol.yold = y; % use y
  sol.fold = f; % use f
end

return
% End of Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
