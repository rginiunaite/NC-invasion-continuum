function tdrTaxis(t, patchId)
%
% function tdrTaxis(t, patchId)
%
% This function computes the taxis approximation on patch patchId for all
% tdr.size equations with approximation (boundary extended) data from
% data{patchId}.y and add to data{patchId}.ydot.
% In this function the taxis limiter function is selected.
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
%*                       data structure before and is removed there)
%*                     - tdrComputeFaceData is called only if
%*                       TDRP.tdr.haveTaxisTerms==true 
%*                     - isTaxis has been renamed isNonZeroVel
%*                     - added non-local term contributions to
%*                       local velocities vx and vy
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

% define limiter functions
limVanLeer = @(theta) (theta+abs(theta)) ./ (1+abs(theta));  % van Leer limiter
limZero    = @(theta) (theta * 0);                           % Zero limiter
limId      = @(theta) (theta);                               % Identity limiter
limiter = limVanLeer;

global TDRP;           % global struct for TDR Problem data

bW = TDRP.grd.boundaryWidth; % get boundary width of patches

if (TDRP.tdr.haveTaxisTerms)
  if (~TDRP.data{patchId}.ComputedFaceData)
    tdrComputeFaceData(t, patchId);
  end
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
    % (1) compute the velocities in the y- or radial direction (vy) 
    vy = zeros(TDRP.grd.ny(patchId)+1, 1);
    isNonZeroVel = false;
    for j=1:TDRP.tdr.size
      if (i==j); % that is done in the diffusion part
	% do nothing
	continue; 
      end;
      switch TDRP.tdr.FTrans.depends(i,j)
       case Zero
        % do nothing
	continue; 
       case Const
	pij = ProbFTrans(i, j, TDRP.params);
	if (pij ~= 0) 
	  vy = vy + pij * TDRP.data{patchId}.uDy(:, j);
	  isNonZeroVel = true;
	end
       case {DependsS, DependsU, DependsSU}
	switch TDRP.tdr.FTrans.depends(i,j)
	 case DependsS
	  pijy = ProbFTrans(i, j, TDRP.params, patchId, 'y');            
	 case {DependsU, DependsSU}
	  pijy = ProbFTrans(i, j, TDRP.params, TDRP.data{patchId}.uAvy, patchId, 'y');      
	 otherwise
	  error('We should never get here!');
	end
	vy = vy + pijy .* TDRP.data{patchId}.uDy(:, j);
	isNonZeroVel = true;
       otherwise
	error('not implemented yet')
      end % switch TDRP.tdr.FTrans.depends(i,j)
    end % for j=1:TDRP.tdr.size
    
    % check if equation i has a non-local term and if yes add the
    % adhesion velocity to vx and vy
    if TDRP.tdr.haveNonLocalTerms
      if (TDRP.tdr.FNonLocal.depends(i) ~= Zero)
	% setup of appropriate matrix G (in cell centres, not boundary).
	switch TDRP.tdr.FNonLocal.depends(i)
	 case {Const}
	  G = ProbFNonLocal(i, TDRP.params)*ones(TDRP.grd.ny(patchId),1);
	 case {DependsS}
	  G = ProbFNonLocal(i, TDRP.params);
	 case {DependsU, DependsSU}
	  G = ProbFNonLocal(i, TDRP.params, ...
			    TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)],:));
	 otherwise
	  error('We should never get here!');
	end
	% eval non-local term using FFT scheme (methodId=2)
	if (size(TDRP.tdr.FNonLocal.R) == [1,1])
	  % single nonlocal term
	  switch TDRP.tdr.FNonLocal.BCs
	   case 'pp'
	    [A] = evalIntegral1D(G, TDRP.tdr.FNonLocal.mask, 2);
	    % add non-local velocity to vy
	    if max(isnan(A(:)))
	      error('brr')
	    end
	    vy = vy + [A(end); A];
	   case 'zz'
	    [A] = evalIntegral1D(G, TDRP.tdr.FNonLocal.mask, 2);
	    if max(isnan(A(:)))
	      error('brr')
	    end
	    % add non-local velocity to vy
	    vy = vy + A;
	    % ***HACK*** enforce zero flux!
	    vy(1)=0;
	    vy(end)=0;
	   case 'vz'
	    LeftLength=TDRP.tdr.FNonLocal.mask.lm+1;
	    Gext = [G(LeftLength:-1:1);G];
	    [A] = evalIntegral1D(Gext, TDRP.tdr.FNonLocal.mask, 2);
	    if max(isnan(A(:)))
	      error('brr')
	    end
	    % add non-local velocity to vy
	    vy = vy + A;
	    % ***HACK*** enforce zero flux!
	    vy(1)=0;
	    vy(end)=0;
	   otherwise
	    error('unsupported value of TDRP.tdr.FNonLocal.BCs');
	  end
	else
	  % single nonlocal term
	  switch TDRP.tdr.FNonLocal.BCs
	   case 'pp'
	    Atotal = zeros(size(G{1}));
	    for jj = 1:size(TDRP.tdr.FNonLocal.R,2)
	      if (TDRP.tdr.FNonLocal.R(i,jj) > 0)
		[A] = evalIntegral1D(G{jj}, TDRP.tdr.FNonLocal.mask{i,jj},2);
		if max(isnan(A(:)))
		  error('brr')
		end
		Atotal = Atotal + A;
	      end
	    end
	    % compute P-function
	    % setup of appropriate matrix P (on grid cell boundaries).
	    switch TDRP.tdr.FNonLocal.depends(i)
	     case {Zero}
	      P = 1*ones(TDRP.grd.ny(patchId)+1,1);
	     case {Const}
	      P = ProbFNonLocal2(i, TDRP.params)*ones(TDRP.grd.ny(patchId)+1,1);
	     case {DependsS}
	      P = ProbFNonLocal(i, TDRP.params);
	     case {DependsU, DependsSU}
	      sol_centre = TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)],:);
	      sol_centre = [
		  sol_centre(end,:)
		  sol_centre
		  sol_centre(1,:)
		  ];
	      sol_bound = 0.5*(sol_centre(1:(end-1),:)+sol_centre(2:end,:));
	      P = ProbFNonLocal2(i, TDRP.params, sol_bound);
	     otherwise
	      error('We should never get here!');
	    end
	    
	    % add non-local velocity to vy
	    vy = vy + P.*[Atotal(end); Atotal];
	   otherwise
	    error('unsupported value of TDRP.tdr.FNonLocal.BCs');
	  end
	end
	isNonZeroVel = true;
      end
    end % of if TDRP.tdr.haveNonLocalTerms
    
    % computation of vy completed
    if (isNonZeroVel)
      % (2) compute taxis/non-local term in y direction
      % for positive velocity
      state = TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], i);
      numer = TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], i) - state;
      denom = state - TDRP.data{patchId}.y(bW-1+[0:TDRP.grd.ny(patchId)], i);
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = (state + 0.5 * (limiter(theta) .* denom)) .* (vy>0) .* (vy);
      % for negative velocity
      state = TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], i);
      % numer is unchangend from above
      denom = TDRP.data{patchId}.y(bW+2+[0:TDRP.grd.ny(patchId)], i) - state;
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = taxisApprox + (state - 0.5 * (limiter(theta) .* denom)) .* (vy<0) .* (vy);
      % add taxisApprox for y direction
      if (TDRP.grd.isAxiSymmetric)
	error('axisymmetric not yet supported in 1D');   
      else
	TDRP.data{patchId}.ydot(:,i) = TDRP.data{patchId}.ydot(:,i) ...
	    - (1.0/TDRP.grd.dy(patchId)) * (taxisApprox(2:end) - taxisApprox(1:end-1));  
      end
    end % if (isNonZeroVel)   
  end % for i=1:TDRP.tdr.size
else % we have a 2D problem
  for i=1:TDRP.tdr.size
    % (1) compute the velocities in the x-direction (vx) and the y
    % or radial direction (vy) 
    vx = zeros(TDRP.grd.ny(patchId)  , TDRP.grd.nx(patchId)+1);
    vy = zeros(TDRP.grd.ny(patchId)+1, TDRP.grd.nx(patchId)  );
    isNonZeroVel = false;
    for j=1:TDRP.tdr.size
      if (i==j); % that is done in the diffusion part
	% do nothing
	continue; 
      end;
      switch TDRP.tdr.FTrans.depends(i,j)
       case Zero
        % do nothing
	continue; 
       case Const
	pij = ProbFTrans(i, j, TDRP.params);
	if (pij ~= 0) 
	  vx = vx + pij * TDRP.data{patchId}.uDx(:, :, j);
	  vy = vy + pij * TDRP.data{patchId}.uDy(:, :, j);
	  isNonZeroVel = true;
	end
       case {DependsS, DependsU, DependsSU}
	switch TDRP.tdr.FTrans.depends(i,j)
	 case DependsS
	  pijx = ProbFTrans(i, j, TDRP.params, patchId, 'x');            
	  pijy = ProbFTrans(i, j, TDRP.params, patchId, 'y');            
	 case {DependsU, DependsSU}
	  pijx = ProbFTrans(i, j, TDRP.params, TDRP.data{patchId}.uAvx, patchId, 'x');      
	  pijy = ProbFTrans(i, j, TDRP.params, TDRP.data{patchId}.uAvy, patchId, 'y');      
	 otherwise
	  error('We should never get here!');
	end
	vx = vx + pijx .* TDRP.data{patchId}.uDx(:, :, j);
	vy = vy + pijy .* TDRP.data{patchId}.uDy(:, :, j);
	isNonZeroVel = true;
       otherwise
	error('not implemented yet')
      end % switch TDRP.tdr.FTrans.depends(i,j)
    end % for j=1:TDRP.tdr.size
    
    % check if equation i has a non-local term and if yes add the
    % adhesion velocity to vx and vy
    if TDRP.tdr.haveNonLocalTerms
      if (TDRP.tdr.FNonLocal.depends(i) ~= Zero)
	% setup of appropriate matrix G (in cell centres, not boundary).
	switch TDRP.tdr.FNonLocal.depends(i)
	 case {Const}
	  G = ProbFNonLocal(i, TDRP.params)*ones(TDRP.grd.ny(patchId),TDRP.grd.nx(patchId));
	 case {DependsS}
	  G = ProbFNonLocal(i, TDRP.params);
	 case {DependsU, DependsSU}
	  G = ProbFNonLocal(i, TDRP.params, ...
			    TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], ...
						 bW+[1:TDRP.grd.nx(patchId)],:));
	 otherwise
	  error('We should never get here!');
	end
	% eval non-local term using FFT2 scheme (methodId=3)
	switch TDRP.tdr.FNonLocal.BCs
	 case 'pp'
	  [A11, A22] = evalIntegral2D_old(G, TDRP.tdr.FNonLocal.mask, 3);
	  % add non-local velocity to vx and vy
	  vx = vx + [A11(:,end) A11];
	  vy = vy + [A22(end,:);A22];
	  if max(isnan(A11(:)))
	    error('brr')
	  end
	  if max(isnan(A22(:)))
	    error('brr')
	  end
	 case {'pppp', 'zzzz'}
	  [A11, A22] = evalIntegral2D(G, TDRP.tdr.FNonLocal.mask, 3);
	  % add non-local velocity to vx and vy
	  vx = vx + A11;
	  vy = vy + A22;
	  if max(isnan(A11(:)))
	    error('brr')
	  end
	  if max(isnan(A22(:)))
	    error('brr')
	  end
	  if strcmp(TDRP.tdr.FNonLocal.BCs, 'zzzz')
	    % enforce zero flux
	    vy(1,:) = 0;
	    vy(end,:) = 0;
	    vx(:,1) = 0;
	    vx(:,end) = 0;
	  end
	 otherwise
	  error('not implemented yet')
	end
	isNonZeroVel = true;
      end
    end % of if TDRP.tdr.haveNonLocalTerms
  
    % computation of vx and vy completed
    if (isNonZeroVel)
      % (2) compute taxis/non-local term in x-direction
      % for positive velocity
      state = TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+  [0:TDRP.grd.nx(patchId)], i);
      numer = TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+1+[0:TDRP.grd.nx(patchId)], i) - state;
      denom = state - TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW-1+[0:TDRP.grd.nx(patchId)], i);
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = (state + 0.5 * (limiter(theta) .* denom)) .* (vx>0) .* (vx);    
      % for negative velocity
      state = TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+1+[0:TDRP.grd.nx(patchId)], i);
      % numer is unchangend from above
      denom = TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)], bW+2+[0:TDRP.grd.nx(patchId)], i) - state;
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = taxisApprox + (state - 0.5 * (limiter(theta) .* denom)) .* (vx<0) .* (vx);
      % add taxisApprox for x direction
      TDRP.data{patchId}.ydot(:,:,i) = TDRP.data{patchId}.ydot(:,:,i) ...
	  - (1.0/TDRP.grd.dx(patchId)) * (taxisApprox(:, 2:end) - taxisApprox(:, 1:end-1));  
      % (3) compute taxis/non-local term in y direction
      % for positive velocity
      state = TDRP.data{patchId}.y(bW+  [0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], i);
      numer = TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], i) - state;
      denom = state - TDRP.data{patchId}.y(bW-1+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], i);
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = (state + 0.5 * (limiter(theta) .* denom)) .* (vy>0) .* (vy);
      % for negative velocity
      state = TDRP.data{patchId}.y(bW+1+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], i);
      % numer is unchangend from above
      denom = TDRP.data{patchId}.y(bW+2+[0:TDRP.grd.ny(patchId)], bW+[1:TDRP.grd.nx(patchId)], i) - state;
      theta = numer ./ (denom - (abs(denom) < 1e-14) ); % avoid division by zero
      taxisApprox = taxisApprox + (state - 0.5 * (limiter(theta) .* denom)) .* (vy<0) .* (vy);
      % add taxisApprox for y direction
      if (TDRP.grd.isAxiSymmetric)
	TDRP.data{patchId}.ydot(:,:,i) = TDRP.data{patchId}.ydot(:,:,i) ...
	    - (1.0/TDRP.grd.dy(patchId)) * (TDRP.data{patchId}.skalYT*taxisApprox(2:end, :) ...
					  - TDRP.data{patchId}.skalYB*taxisApprox(1:end-1, :));   
      else
	TDRP.data{patchId}.ydot(:,:,i) = TDRP.data{patchId}.ydot(:,:,i) ...
	    - (1.0/TDRP.grd.dy(patchId)) * (taxisApprox(2:end, :) - taxisApprox(1:end-1, :));  
      end
    end % if (isNonZeroVel)
  end % for i=1:TDRP.tdr.size
end % if TDRP.grd.is1D

return
% end of function
