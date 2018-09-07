function pij = ProbFTrans(i, j, params, varargin)
% function pij = ProbFTrans(i, j, params, varargin)
%
% This function is part of the DiffReac TDR model.
%
% This function implements the diffusion coefficient and taxis
% coefficient functions.
% Use ComparisonAxisymmetricAnd2D() to run the simulation.
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/ProbFTrans.m
%* Date created  : 2006, January 25
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

if ( (i == j) & (i == 1))
  pij = params.diffConst(1);
else  
  error('ProbFTrans:: unimplemented.');
end

return
% end of function
