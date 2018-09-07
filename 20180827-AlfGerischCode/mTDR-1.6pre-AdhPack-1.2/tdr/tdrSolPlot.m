function [ymin, ymax]=tdrSolPlot(y, EqNo)
%
% function [ymin, ymax]=tdrSolPlot(y, EqNo)
%
% This function plots the solution component EqNo of the solution of the 
% TDR system provided in vector y (in standard order) in the current axis. 
%
% Returned are the minimum and the maximum value of the solution of 
% equation EqNo.
%
% This function calls: tdrUtil/plotPatch()
%
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdr/tdrSolPlot.m
%* Date created  : 2005, September 25
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

				
global TDRP;           % global struct for TDR Problem data
				 
ymin=min(y(EqNo:TDRP.tdr.size:end)); % compute minimum solution value
ymax=max(y(EqNo:TDRP.tdr.size:end)); % compute maximum solution value 


if ~TDRP.grd.is1D
  % plotparams is a struct to control plotPatch()
  %plotparams.methodId = 1;
  plotparams.methodId = 2;

  for patchId = 1:TDRP.grd.nop
    % Determine ystart and yend such that y(ystart:yend) corresponds to the
    % current patchId
    ystart = TDRP.grd.ps(patchId)-1+EqNo;
    yend   = TDRP.grd.pe(patchId)-1+EqNo;
    % extract solution values from y to zvals
    zvals = reshape(y(ystart:TDRP.tdr.size:yend), ...
		    TDRP.grd.ny(patchId), TDRP.grd.nx(patchId));
    % plot the current patch
    plotPatch(TDRP.grd.x0(patchId),TDRP.grd.dx(patchId),TDRP.grd.nx(patchId), ...
	      TDRP.grd.y0(patchId),TDRP.grd.dy(patchId),TDRP.grd.ny(patchId), ...
	      zvals, plotparams); 
    hold('on')
  end
else
  plot(TDRP.grd.cellCentreMatrix, y(EqNo:TDRP.tdr.size:end))
end

return;
% end of function
