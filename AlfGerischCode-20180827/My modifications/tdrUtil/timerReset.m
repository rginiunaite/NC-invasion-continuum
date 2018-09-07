function timerReset()
%
% function tdrReset()
%
% This function resets all timer and counter in global struct TIMER to zero.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrUtil/timerReset.m
%* Date created  : 2005, September 29
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     :
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2005 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************

global TIMER;          % global struct for timer

% reset all timer
TIMER.tdrReacV  = 0.0;
TIMER.tdrReac   = 0.0;
TIMER.tdrPrepTrans = 0.0;
TIMER.tdrFinTrans = 0.0;
TIMER.tdrDiff   = 0.0;
TIMER.tdrTaxis  = 0.0;

% reset all counter
TIMER.tdrReacVC  = 0;
TIMER.tdrReacC   = 0;
TIMER.tdrPrepTransC =0;
TIMER.tdrFinTransC = 0;
TIMER.tdrDiffC   = 0;
TIMER.tdrTaxisC  = 0;

return
% end of function