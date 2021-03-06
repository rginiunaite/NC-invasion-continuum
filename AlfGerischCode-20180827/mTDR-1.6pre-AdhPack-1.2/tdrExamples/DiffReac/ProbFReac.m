function [ydot] = ProbFReac(t, coords, y, params)
% function ydot = ProbFReac(t, coords, y, params)
%
% This function is part of the DiffReac TDR model.
%
% This function implements the reaction function in a vectorized manner.
% The arguments for the reaction terms are (coords(1,i), coords(2,i),
% y(:,i)) for i=1:size(y,2). Return value ydot has the same size as y.
% Use ComparisonAxisymmetricAnd2D() to run the simulation.
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/ProbFReac.m
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

if strcmp(params.selectReaction, 'none')
  ydot = 0*y;
elseif strcmp(params.selectReaction, 'logistic')
  alpha = params.ReactionParams(1);
  beta  = params.ReactionParams(2);
  ydot  = alpha * ((y.^beta).*(1-y));
elseif strcmp(params.selectReaction, 'genlogistic')
  alpha = params.ReactionParams(1);
  beta  = params.ReactionParams(2);
  ydot = -alpha * (y.*(1-y).*(beta-y));
else
  error('ProbFReac: unknown value of params.selectReaction');
end

return;