function tdrDiff(t, patchId)
%
% function tdrDiff(t, patchId)
%
% This function computes the diffusion approximation on patch patchId for all 
% tdr.size equations with approximation (boundary extended) data from 
% data{patchId}.y and add to data{patchId}.ydot.
%
% This function calls external functions:
%        tdrComputeFaceData(), ProbFTrans()
% This function calls internal functions:
%        none
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdr/tdrDiff.m
%* Date created  : 2005, September 25
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.1 (AG, 2006 06 08)
%*                     - derivatives in the radial direction are now
%*                       computed correctly in the axisymmetric case.
%*                 1.2 (AG, 2007 01 03)
%*                     - added patchId and 'x' or 'y' to the end of the list  
%*                       of arguments, when calling ProbFTrans.m (unless
%*                       the transport function is a constant). This
%*                       enables space dependent parameters in the
%*                       transport functions.
%*                 1.3 (AG, 2008 09 08)
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

if (~TDRP.data{patchId}.ComputedFaceData)
  tdrComputeFaceData(t, patchId);
end
                
%%% CODES %%%%
  Zero       = 0;                              % f \equiv 0
  Const      = 1;                              % f \equiv const
  DependsT   = 2;                              % f = f(t)
  DependsS   = 4;                              % f = f(x,y)
  DependsU   = 8;                              % f = f(u)
  DependsTS  = DependsT + DependsS;            % f = f(t, x, y)
  DependsTU  = DependsT + DependsU;            % f = f(t, u)
  DependsSU  = DependsS + DependsU;            % f = f(x, y, u)
  DependsTSU = DependsT + DependsS + DependsU; % f = f(t, x, y, u)

if TDRP.grd.is1D
  for i=1:TDRP.tdr.size
    switch TDRP.tdr.FTrans.depends(i,i)
     case  Zero
      % do nothing
      continue
     case Const
      pii = ProbFTrans(i, i, TDRP.params);
      % y-diffusion (that is radial diffusion (not supported yet))
      % for equation i 
      ydiff = (pii/TDRP.grd.dy(patchId)) * ...
	      (TDRP.data{patchId}.uDy(2:end, i) - TDRP.data{patchId}.uDy(1:end-1, i));
     case {DependsS, DependsU, DependsSU}
      switch TDRP.tdr.FTrans.depends(i,i)
       case DependsS
	pii = ProbFTrans(i, i, TDRP.params, patchId, 'y');            
       case {DependsU, DependsSU}
	pii = ProbFTrans(i, i, TDRP.params, TDRP.data{patchId}.uAvy, patchId, 'y');      
       otherwise
	error('We should never get here!');
      end
       % y-diffusion (that is radial diffusion (not supported yet)) for equation i
      ydiff = (1.0/TDRP.grd.dy(patchId)) * ...
	      (  pii(2:end  ) .* TDRP.data{patchId}.uDy(2:end  , i) ...
	       - pii(1:end-1) .* TDRP.data{patchId}.uDy(1:end-1, i));
     otherwise
      error('not implemented yet')
    end % switch TDRP.tdr.FTrans.depends(i,i)
    % add to TDRP.data{patchId}.ydot
    TDRP.data{patchId}.ydot(:,i) = ...
	TDRP.data{patchId}.ydot(:,i) + ydiff;
  end % for i=1:TDRP.tdr.size
else
  for i=1:TDRP.tdr.size
    switch TDRP.tdr.FTrans.depends(i,i)
     case  Zero
      % do nothing
      continue
     case Const
      pii = ProbFTrans(i, i, TDRP.params);
      % x-diffusion for equation i
      xdiff = (pii/TDRP.grd.dx(patchId)) * ...
	      (TDRP.data{patchId}.uDx(:, 2:end, i) - TDRP.data{patchId}.uDx(:, 1:end-1, i));
      % y-diffusion (that is radial diffusion) for equation i
      if (TDRP.grd.isAxiSymmetric)
	ydiff = (pii/TDRP.grd.dy(patchId)) * ...
		(TDRP.data{patchId}.skalYT*TDRP.data{patchId}.uDy(2:end, :, i) ...
		-TDRP.data{patchId}.skalYB*TDRP.data{patchId}.uDy(1:end-1, :, i));
      else
	ydiff = (pii/TDRP.grd.dy(patchId)) * ...
		(TDRP.data{patchId}.uDy(2:end, :, i) - TDRP.data{patchId}.uDy(1:end-1, :, i));
      end
     case {DependsS, DependsU, DependsSU}
      % x-diffusion for equation i
      switch TDRP.tdr.FTrans.depends(i,i)
       case DependsS
	pii = ProbFTrans(i, i, TDRP.params, patchId, 'x');            
       case {DependsU, DependsSU}
	pii = ProbFTrans(i, i, TDRP.params, TDRP.data{patchId}.uAvx, patchId, 'x');      
       otherwise
	error('We should never get here!');
      end
      xdiff = (1.0/TDRP.grd.dx(patchId)) * ...
	      (  pii(:, 2:end  ) .* TDRP.data{patchId}.uDx(:, 2:end  , i) ...
	       - pii(:, 1:end-1) .* TDRP.data{patchId}.uDx(:, 1:end-1, i));
      % y-diffusion (that is radial diffusion) for equation i
      switch TDRP.tdr.FTrans.depends(i,i)
       case DependsS
	pii = ProbFTrans(i, i, TDRP.params, patchId, 'y');            
       case {DependsU, DependsSU}
	pii = ProbFTrans(i, i, TDRP.params, TDRP.data{patchId}.uAvy, patchId, 'y');      
       otherwise
	error('We should never get here!');
      end
      if (TDRP.grd.isAxiSymmetric)
	ydiff = (1.0/TDRP.grd.dy(patchId)) * ...
		(  pii(2:end  , :) .* (TDRP.data{patchId}.skalYT*TDRP.data{patchId}.uDy(2:end  , :, i)) ...
		 - pii(1:end-1, :) .* (TDRP.data{patchId}.skalYB*TDRP.data{patchId}.uDy(1:end-1, :, i)));
      else
	ydiff = (1.0/TDRP.grd.dy(patchId)) * ...
		(  pii(2:end  , :) .* TDRP.data{patchId}.uDy(2:end  , :, i) ...
		 - pii(1:end-1, :) .* TDRP.data{patchId}.uDy(1:end-1, :, i));
      end
     otherwise
      error('not implemented yet')
    end % switch TDRP.tdr.FTrans.depends(i,i)
    % add to TDRP.data{patchId}.ydot
    TDRP.data{patchId}.ydot(:,:,i) = ...
	TDRP.data{patchId}.ydot(:,:,i) + xdiff + ydiff;
  end % for i=1:TDRP.tdr.size
end % if TDRP.grd.is1D
  
return
% end of function
