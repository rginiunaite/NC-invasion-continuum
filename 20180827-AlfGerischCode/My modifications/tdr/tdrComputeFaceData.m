function tdrComputeFaceData(t, patchId)
%
% function tdrComputeFaceData(t, patchId)
%
% This function is called for each patch and computes averaged data and 
% derivatives on the grid faces of each cell in that patch if 
% TDRP.data{patchId}.ComputedFaceData has the value false; otherwise nothing is 
% done. The computed data is stored in TDRP.data{patchId} fields as described 
% below and TDRP.data{patchId}.ComputedFaceData is set to true.
% The following four fields are created. Each is a 3D array with third 
% dimension equal to tdr.size and the first two dimensions according to 
% the number of horizontal/vertical grid faces in the patch. 
%   data{patchId}.uDx    ...  x-derivative approximation of all solution components 
%                    on the right cell boundary of each cell for patch 
% 		     patchId (incl. on the left boundary of the patch) 
%   data{patchId}.uDy    ...  y-derivative approximation of all solution components 
%                    on the upper cell boundary of each cell for patch 
% 		     patchId (incl. on the bottom boundary of the patch)
% 		     Note: data.uDy is not scaled in order to account for
%                    axial symmetry if TDRP.grd.isAxiSymmetric == true.
%   data{patchId}.skalYT ...  scaling matrix for derivative in radial (y) direction
%                    to be applied on top cell faces of current patch
%                    (only created if TDRP.grd.isAxiSymmetric == true).
%   data{patchId}.skalYB ...  scaling matrix for derivative in radial (y) direction
%                    to be applied on bottom cell faces of current patch
%                    (only created if TDRP.grd.isAxiSymmetric == true).
%   data{patchId}.uAvx   ...  approximations to all solution components on the right 
%                    cell boundary of each cell for patch patchId (incl. 
% 		     on the left boundary of the patch) by averaging.
%   data{patchId}.uAvy   ...  approximations to all solution components on the upper 
%                    cell boundary of each cell for patch patchId (incl. 
%      		     on the bottom boundary of the patch) by averaging.
%
% This function calls external functions:
%        none
% This function calls internal functions:
%        none
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdr/tdrComputeFaceData.m
%* Date created  : 2005, September 25
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.1 (AG, 2006 06 08)
%*                     - added data.skalYT and data.skalYB to correctly
%*                       compute derivatives in radial direction.
%*                 1.2 (AG, 2008 09 08)
%*                     - TDRP.grd.boundaryWidth is a new field of
%*                       the grd structure (has been part of the
%*                       data structure before and is removed there 
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2008 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************


global TDRP;           % global struct for TDR Problem data
bW = TDRP.grd.boundaryWidth; % get boundary width of patches

if (t ~= TDRP.data{patchId}.t)
  error('tdrComputeFaceData::inconsistent time points.')
end

if TDRP.data{patchId}.ComputedFaceData
  return
end

if (TDRP.grd.is1D)
  if (patchId ~= 1)
    error('1D only for one patch.');
  end
  % Compute derivative approximation of solution on cell boundaries
  % and store in TDRP.data{1}.uDy 
    TDRP.data{patchId}.uDy ...
      = (1/TDRP.grd.dy(patchId)) * ( ...
	  TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], :) ...
	 -TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], :) );
  % Compute averages of concentration in cell faces in TDRP.data{patchId}.uAvy
  TDRP.data{patchId}.uAvy = ...
      (1/2) * ( ...
	  TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], :) ...
	 +TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], :) );
  % set computed flag to true
  TDRP.data{patchId}.ComputedFaceData = true;
else
  % Compute x and y derivative approximation of solution on cell 
  % boundaries for patch patchId and store in TDRP.data{patchId}.uDx and TDRP.data{patchId}.uDy
  TDRP.data{patchId}.uDx ...
      = (1/TDRP.grd.dx(patchId)) * ( ...
	  TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+1+[0:TDRP.grd.nx(patchId)], :) ...
	 -TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+  [0:TDRP.grd.nx(patchId)], :) );
  TDRP.data{patchId}.uDy ...
      = (1/TDRP.grd.dy(patchId)) * ( ...
	  TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], :) ...
	 -TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], :) );
  if (TDRP.grd.isAxiSymmetric)
    % TODO: THIS IS TIME AND SOL INDEPENDENT STUFF AND SO SHOULD BE
    % MOVED ELSEWHERE!!
    % compute cell centers in y-direction (including the centre of the cell
    % below the current patch)
    cellCentresR = TDRP.grd.y0(patchId) + ([0:TDRP.grd.ny(patchId)]-0.5)* TDRP.grd.dy(patchId);
    mdim = size(TDRP.data.uDy,1)-1;
    % compute scaling matrices for derivative in radial (y) direction on
    % top (skalYT) and bottom (skalYB) cell faces 
    TDRP.data{patchId}.skalYT = spdiags(1+(TDRP.grd.dy(patchId)/2)./cellCentresR(2:end)', 0, mdim, mdim);
    TDRP.data{patchId}.skalYB = spdiags(1-(TDRP.grd.dy(patchId)/2)./cellCentresR(1:end-1)', 0, mdim, mdim);
  end
  % Compute averages of concentration in x-cell faces in TDRP.data{patchId}.uAvx 
  % and on y-cell faces in TDRP.data{patchId}.uAvy
  TDRP.data{patchId}.uAvx ...
      = (1/2) * ( ...
	  TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+1+[0:TDRP.grd.nx(patchId)], :) ...
	 +TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+  [0:TDRP.grd.nx(patchId)], :) );
  TDRP.data{patchId}.uAvy = ...
      (1/2) * ( ...
	  TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], :) ...
	 +TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], :) );
  % set computed flag to true
  TDRP.data{patchId}.ComputedFaceData = true;
end

return;

% end of function
