function [BCtype,BCval] = ProbBCs(patchId, bdry, t, xvals, yvals, params)
% function [BCtype,BCval] = ProbBCs(patchId, bdry, t, xvals, yvals, params)
%
% This function is part of the linear advection model problem, 
% which illustrates the basic usage of the Matlab TDR system.
%
% This function implements the boundary conditions. Use
% LinearAdvection() to run the simulation.
% 
% The full model consists of the five m-files: LinearAdvection.m, 
% ProbBCs.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/LinearAdvection/ProbBCs.m
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of Standard Codes               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defined boundaries
left   = 1;
right  = 2;
bottom = 3;
top    = 4;
% defined boundary condition types
None      = 0;
ZeroFlux  = 1;
Dirichlet = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of Definition of Standard Codes        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This example has two equations. For the first equation 
% (containing the advection term) we prescribe Dirichlet 
% zero BCs along the inflow edge and zero flux elsewhere,
% for the second (provides the constant velocity profile 
% only) no BCs are prescribed mathematically but we need 
% derivatives on the boundary so simply use the  
% appropriate Dirichlet BC.

if strcmp(params.direction,'x')
  variableDir = xvals;
  if (params.velocity > 0)
    inflow = left;   % left boundary is inflow boundary
  else
    inflow = right;  % right boundary is inflow boundary
  end
else
  variableDir = yvals;
  if (params.velocity > 0)
    inflow = bottom; % bottom boundary is inflow boundary
  else
    inflow = top;    % top boundary is inflow boundary
  end
end

vec1 = ones(size(xvals));
if (bdry == inflow)
  BCtype = vec1 * [Dirichlet Dirichlet];
  BCval  = [vec1*0  params.velocity*variableDir];
else
  BCtype = vec1 * [ZeroFlux Dirichlet];
  BCval  = [vec1*NaN  params.velocity*variableDir];
end

return;
% end of function
