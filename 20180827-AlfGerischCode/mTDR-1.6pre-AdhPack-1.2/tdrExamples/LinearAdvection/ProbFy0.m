function y0 = ProbFy0(coords, params)
% function y0 = ProbFy0(coords, params)
%
% This function is part of the linear advection model problem, 
% which illustrates the basic usage of the Matlab TDR system.
%
% This function implements the initial conditions. 
% Use LinearAdvection() to run the simulation.
% 
% The full model consists of the five m-files: LinearAdvection.m, 
% ProbBCs.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/LinearAdvection/ProbFy0.m
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

% Return the initial values of all species in (coords(1), coords(2)).
if strcmp(params.direction,'x')                           
 variableDir = coords(1);
else
 variableDir = coords(2);
end

% y0(1) is a profile constant in the direction perpendicular to 
% params.direction
y0(1) = max(0, 1 - 10*abs(0.5-variableDir));


% y0(2) has a constant slope of params.velocity in the direction 
% params.direction
y0(2) = params.velocity * variableDir;

return;
% end of function