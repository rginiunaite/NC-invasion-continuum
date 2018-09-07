function pij = ProbFTrans(i, j, params, varargin)
% function pij = ProbFTrans(i, j, params, varargin)
%
% This function is part of the linear advection model problem, 
% which illustrates the basic usage of the Matlab TDR system.
%
% This function implements the diffusion coefficient and taxis
% coefficient functions. Use LinearAdvection() to run the simulation.
% 
% The full model consists of the five m-files: LinearAdvection.m, 
% ProbBCs.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/LinearAdvection/ProbFTrans.m
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


% Return function values of parameter function p_{i,j}(...)
if (i == j)
    error('unimplemented');
else 
  if ((i==1) && (j==2))
    pij = 1;
  else 
    error('unimplemented');
  end
end

return
% end of function
