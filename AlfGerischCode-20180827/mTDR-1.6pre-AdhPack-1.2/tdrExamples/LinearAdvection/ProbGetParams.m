function [grd, tdr, params] = ProbGetParams(velocity, direction)
% function [grd, tdr, params] = ProbGetParams(velocity, direction)
%
% This function is part of the linear advection model problem, 
% which illustrates the basic usage of the Matlab TDR system.
%
% This function sets up the parameters of the model and those for the
% simulation. Parameter direction is either 'x' or 'y' and prescribes 
% the direction of transport. Parameter velocity is a real number and
% gives the (positive or negative) speed of propagation.
% Use LinearAdvection() to run the simulation.
% 
% The full model consists of the five m-files: LinearAdvection.m, 
% ProbBCs.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/LinearAdvection/ProbGetParams.m
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

%
% initialise struct params
%
params.velocity = velocity;   % negative for move to the left (bottom) edge
params.direction= direction;  % set direction of movement
%
% params struct completed
%

%
% initialise struct tdr
%
tdr.size = 2;          % number of species (PDEs)
tdr.tvec = 0:0.1:1;    % vector of output times      

%%% CODES %%%%
  Zero       = 0;                              % f \equiv 0
  Const      = 1;                              % f \equiv const
  DependsT   = 2;                              % f = f(t)
  DependsS   = 4;                              % f = f(x,y)
  DependsU   = 8;                              % f = f(u)
  DependsTS  = DependsT + DependsS;            % f = f(t, x, y)
  DependsTU  = DependsT + DependsU;            % f = f(t, u)
  DependsSU  = DependsS + DependsU;            % f = f(x, y, u)
  DependsTSU = DependsT + DependsS + DependsU; % f = f(t, x, y, u)

% describe reaction function ProbFReac.m 
% if only Zero entries then ProbFReac.m is never called
tdr.FReac.depends = [ Zero; Zero ];  
tdr.FReac.vec = true;  % true if ProbFReac.m can be called vectorized

% describe transport function ProbFTrans.m
% if only Zero entries then ProbFTrans.m is never called
tdr.FTrans.depends = [Zero    Const    
                      Zero    Zero]; 
%
% tdr struct completed
%

%
% define the spatial domain via its patches and store data in grd struct
%
%%%ONE BIG PATCH (unit square)%%%
if (1)
  grd.isAxiSymmetric = false;
  grd.nop = 1;             % Number Of Patches
  grd.ngb = [0 0 0 0];     % patchId of [left right bottom top] neighbour patch
			   % 0 means no neighbour in that direction
  if strcmp(params.direction,'x')                           
    grd.nx = [100];          % number of cells in x-direction in patch
    grd.ny = [10];           % number of cells in y-direction in patch
  elseif strcmp(params.direction,'y')
    grd.nx = [10];           % number of cells in x-direction in patch
    grd.ny = [100];          % number of cells in y-direction in patch
  else
    error('Provided value for direction not allowed.');
  end
  grd.dx = 1/grd.nx;       % x-cell width in patch
  grd.dy = 1/grd.ny;       % y-cell width in patch
  grd.x0 = 0;              % x-coordinate of lower left corner of patch
  grd.y0 = 0;              % y-coordinate of lower left corner of patch
end
%%%End of ONE BIG PATCH%%%

%
% grd struct completed
%

return;
% end of function
