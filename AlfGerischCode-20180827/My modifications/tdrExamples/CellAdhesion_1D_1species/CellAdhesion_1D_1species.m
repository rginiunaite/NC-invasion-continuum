function [t,y]=CellAdhesion_1D_1species(varargin)
% function [t,y]=CellAdhesion_1D_1species(varargin)
%
% This is the main function to run a simulation and can just be
% called with one argument (default values given below):
%    params.domainlength    = 1;
%    params.gridCells       = 400;
%    params.BCs             = 'vz'; % or 'pp'
%    inputStruct.params = params;
%    CellAdhesion_1D_1species(inputStruct);
% 
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : 
%* Date created  : 2008, Sep 30
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.0 initial version (Alf Gerisch)
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2006-2008 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%******************************************************************************


clear functions
format compact;

oldpath=path();          % save current Matlab path
% include mTDR base directory in path and setup mTDR system
path('../../../mTDR-1.6pre-AdhPack-1.2', path); 
tdrSetup;
% add rowmap to Matlab path
path('~/Matlab/rowmap-0.93',path); 


% print TDR, rowmap and problem directories used 
% as convenience for the user
which tdrInit            % is in tdr/tdr directory
which rowmap
which ProbGetParams      % is in problem directory

% Select time integration scheme
%integrator = 'rowmap';
integrator = 'ode45';
% Set tolerance for time integration
tol=1e-7;
% Set output function
clear function localOutputfun; % to clear persistent variables
outputfun=@localOutputfun; % This function is defined below.

% init the TDR problem (initialises the data in TDRP and returns 
% the initial value y0 and the vector of output times tspan)
clear function ProbFy0; % to reinitialize persistent variables
[y0, tspan] = tdrInit(varargin{:});

timerReset();           % reset all timer
switch integrator
 case 'ode45'
  disp('Running ode45 ...');
  % install output function for odesolver
  clear options;
  options=[];
  options = odeset(options, 'OutputFcn', outputfun);
  options = odeset(options, 'AbsTol', tol, 'RelTol', tol);
  tic
  [t,y] = ode45(@tdrFdgl, tspan, y0, options);
  toc
 case 'rowmap'
  disp('Running rowmap ...');
  clear options;
  options.OptWarnMiss = 1;
  options.AbsTol = tol;
  options.RelTol = tol;
  options.RHSnonautonomous = 0; % the MOL-ODE system is explicitely
				% NOT depending on t. 
  options.OutputFun = outputfun;
  options.MaxKrylovSteps = 50;
  options.MaxNumberOfSteps = 5000000;
  %options.IncludeGridPointsInDenseOutput = 1;
  options.InitialStep = 1e-5;
  options.ReturnMode = 2;
  options.OutputCallMode = 2;
  options.ContOutputInterp = 1;
  options
  tic
  [t,y,stat, RowmapNextInitialStep] = ...
	  rowmap(@tdrFdgl, tspan, y0, options);
  toc
  rowmapidid = stat(1);
  if (rowmapidid ~= 1)
    timerPrint('main::');   % print current values of timer
    error(['ROWMAP terminated with IDID = ' num2str(rowmapidid)]);
  end
 otherwise
  error('unknown integrator');
end
%%%%%%%%%%%%%%%%%%%%%%%%
timerPrint('main::');   % print current values of timer


%{
global TDRP;           % global struct for TDR Problem data
figure(2)
clf;
set(gcf,'DoubleBuffer','on');
grid.x2 = TDRP.grd.yc{1};
waterfall(grid.x2,t,y);
view([25 60])
mm=[min(y') 
    max(y')];
format short e
disp(mm')
%}

% reset Matlab path
path(oldpath);
  
return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = localOutputfun(t, y, flag)
global TDRP;           % global struct for TDR Problem data

persistent saved;
saveEvery = 25;
doPlot=false;
doPlot=true;

res = 0;

if (strcmp(flag,'init'))
  disp('Received flag ''init'' in output function.');
  saved.params = TDRP.params;
  saved.tdr    = TDRP.tdr;
  saved.grd    = TDRP.grd;
  saved.t      = t(1);
  
  if (size(y,2)~=[1])
    error('wrong size of y')
  end
  saved.y = y';
   
  save('saved.mat', 'saved');
  disp(['Saved to saved.mat at time t = ' num2str(t(1))]);
  if (doPlot)
    plotfun(saved.t(end), saved.y(end,:));
    %disp('Hit button to continue...');
  end
  
  return
end

if (strcmp(flag,'done'))
  disp('Received flag ''done'' in output function.');
  save('saved.mat', 'saved');
  disp(['saved.mat .']);
  return
end

disp(['Reached time t = ' num2str(t)]);

if (size(y,2)~=[1])
  error('wrong size of y')
end
if (size(t) ~= [1 1])
  error('wrong size of t');
end

itmp=length(saved.t)+1;
saved.t(itmp) = t;
saved.y(itmp,:) = y';

if (~mod(length(saved.t), saveEvery))
  save('saved.mat', 'saved');
  disp(['Saved to saved.mat at time t = ' num2str(t(end))]);
end

if (doPlot)
    plotfun(saved.t(end), saved.y(end,:));
end

return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = plotfun(t, y)

global TDRP

%disp([t min(y) max(y), sum(y)]);
figure(1)
clf;
%set(gcf, 'Position', [0 0 600 800]);
%set(gcf,'renderer','zbuffer');
set(gcf,'DoubleBuffer','on');

%grid.x1 = TDRP.grd.xc{1};
grid.x2 = TDRP.grd.yc{1};

% solution plot
plot(grid.x2,y);
axis([0 1 -0.1 1.1])
drawnow
pause(0.01)

res = 0;
return;
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

