function timerPrint(idstr)
%
% function tdrReset(idstr)
%
% This function prints the current value of all timer and counter in the  
% global struct TIMER. Each printed line starts with the string 'idstr'
% provided as argument to the function.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrUtil/timerPrint.m
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

disp([idstr 'Timer'])
disp([idstr 'TIMER.tdrReacV     = ' num2str(TIMER.tdrReacV)     ...
      ' (' num2str(TIMER.tdrReacVC)     ' calls)']);
disp([idstr 'TIMER.tdrReac      = ' num2str(TIMER.tdrReac)      ...
      ' (' num2str(TIMER.tdrReacC)      ' calls)']);
disp([idstr 'TIMER.tdrPrepTrans = ' num2str(TIMER.tdrPrepTrans) ...
      ' (' num2str(TIMER.tdrPrepTransC) ' calls)']);
disp([idstr 'TIMER.tdrDiff      = ' num2str(TIMER.tdrDiff)      ...
      ' (' num2str(TIMER.tdrDiffC)      ' calls)']);
disp([idstr 'TIMER.tdrTaxis     = ' num2str(TIMER.tdrTaxis)     ...
      ' (' num2str(TIMER.tdrTaxisC)     ' calls)']);
disp([idstr 'TIMER.tdrFinTrans  = ' num2str(TIMER.tdrFinTrans)  ...
      ' (' num2str(TIMER.tdrFinTransC)  ' calls)']);

return
% end of function
