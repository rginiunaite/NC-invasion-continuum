function [BCtype,BCval] = ProbBCs(patchId, boundary, t, xvals, yvals, params)
% function [BCtype,BCval] = ProbBCs(patchId, boundary, t, xvals, yvals, params)
%
% This function is part of the DiffReac TDR model.
%
% This function implements the boundary conditions.
% Use ComparisonAxisymmetricAnd2D() to run the simulation.
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/ProbBCs.m
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of Standard Codes               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Codes to identify the boundaries of a cell or patch:
  left       = 1;
  right	     = 2;
  bottom     = 3;
  top	     = 4;

% Codes to define boundary condition types
  None       = 0;
  ZeroFlux   = 1;
  Dirichlet  = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of Definition of Standard Codes        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This example has one equation and mainly zero flux bounadry conditions:
BCtype = ZeroFlux * ones(size(xvals));
BCval  = NaN      * ones(size(xvals));

if (params.selectICandBC == 1)
  if (boundary == top)
    BCtype = Dirichlet * ones(size(xvals));
    BCval  = 1.0       * ones(size(xvals));
  end
elseif (params.selectICandBC == 2)
  % do nothing
else
  error('ProbBCs: unknown value of params.selectICandBC');
end

return;
% end of function
