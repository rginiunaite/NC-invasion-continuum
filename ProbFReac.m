function [ydot] = ProbFReac(t, coords, y, params)
% function ydot = ProbFReac(t, coords, y, params)
%
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

[Lt,Ltdot] = domaingrowth(t,params);

ydot = [
    -Ltdot/Lt * y(params.eq.u,:)
    -Ltdot/Lt * y(params.eq.v,:)
    -Ltdot/Lt * y(params.eq.c,:) ...
      - params.reac_degradation_alpha*(y(params.eq.u,:)+y(params.eq.v,:)).*y(params.eq.c,:) ...
      + params.reac_production_kappa * y(params.eq.c,:) .* (1-y(params.eq.c,:))
];

return;