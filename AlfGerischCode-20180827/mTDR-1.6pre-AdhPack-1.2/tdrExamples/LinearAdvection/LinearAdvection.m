function [t,y] = LinearAdvection()
% function [t,y] = LinearAdvection()
%
% This function is part of the linear advection model problem, 
% which illustrates the basic usage of the Matlab TDR system.
%
% This is the main function to run a simulation and can just be 
% called without arguments.
% 
% The full model consists of the five m-files: LinearAdvection.m, 
% ProbBCs.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/LinearAdvection/LinearAdvection.m
%* Date created  : 2006, January 27
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     :
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2006 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************

clear functions
format compact;

oldpath=path();          % save current Matlab path
% include ./tdr and ./tdrUtil in path
path('../../tdrUtil',path); 
path('../../tdr',path); 

% print TDR directory and problem directory used 
% as convenience for the user
which tdrInit            % is in tdr directory
which ProbGetParams      % is in problem directory
    
timerReset();           % reset all timer

% set param and options for ODE solver as empty
param = [];
options = [];
% install output function for odesolver (see tdrUtil subdir)
options = odeset('OutputFcn', @outputfunPlot);

% init the TDR problem (initialises the data in TDRP and returns 
% the initial value y0 and the vector of output times tspan)
velocity = -0.2
direction = 'y'
[y0, tspan] = tdrInit(velocity, direction);

% plot initial conditions
outputfunPlot(tspan(1), y0, '', []);
%disp('Press any key to continue ...'); pause
disp('Running...');

% call to ODE solver performs the time integration
[t,y] = ode45(@tdrFdgl, tspan, y0, options, param);

timerPrint('main::');   % print current values of timer

disp('Press any key to continue ...'); pause

% show simulation results
for ii=1:length(t)
  clf;
  outputfunPlot(t(ii), y(ii,:)', '', []);
  view(2)
  pause(0.1)
end

% reset Matlab path
path(oldpath);
  
return
