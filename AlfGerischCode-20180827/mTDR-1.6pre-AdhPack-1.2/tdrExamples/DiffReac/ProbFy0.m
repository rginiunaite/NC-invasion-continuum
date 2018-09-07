function y0 = ProbFy0(coords, params)
% function y0 = ProbFy0(coords, params)
%
% This function is part of theDiffReac TDR model.
%
% This function implements the initial conditions.
% Use ComparisonAxisymmetricAnd2D() to run the simulation.
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/ProbFy0.m
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

if (params.selectICandBC == 1)
  y0(1) = (coords(2)-params.yshift)^2 * ...
	  exp( -50*((coords(2)-params.yshift)-1.0).^2 );
elseif (params.selectICandBC == 2)
  y0(1) = exp( -50*((coords(2)-params.yshift)-0.3).^2);
else
  error('ProbFy0: unknown value of params.selectICandBC');
end

return
% end of function
