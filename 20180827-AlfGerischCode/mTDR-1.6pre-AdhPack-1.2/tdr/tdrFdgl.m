function  ydot = tdrFdgl(t, y, param)
%
% function ydot = tdrFdgl(t, y, param)
%
% This is the main function which provides the right-hand side F(t,y) of 
% the MOL-ODE of the TDR system. Here t is the current time and y is the 
% current solution approximation of the MOL-ODE system (in standard order).
% Parameter param is the struct variable param as provided in the call 
% to the ODE solver and not the struct params created in ProbGetParams().
% 
% The return value ydot is F(t,y) in standard order.
% 
% This function is also used for timing using the global struct TIMER.
%  
% This function calls external functions: 
%        tdrDiff(), tdrTaxis(), ProbFReac()
% This function calls internal functions: 
%        tdrPrepTrans()
% 
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdr/tdrFdgl.m
%* Date created  : 2005, September 25
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.1 (AG, 2008 09 08)
%*                     - TDRP.grd.boundaryWidth is a new field of
%*                       the grd structure (has been part of the
%*                       TDRP.data structure before and is removed
%*                       there). It is not defined here anymore.  
%*                     - Make use of TDRP.tdr.have{Reaction,Diffusion,
%*                       Taxis,NonLocal}Terms here and do not
%*                       define your own version of those. 
%*                     - Subfunction tdrPrepTrans() does not need
%*                       bW argument anymore, i.e. removed. 
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
global TIMER;          % global struct for timer

% init ydot to zero and size of y
ydot = zeros(size(y));


% loop over all patches
for patchId = 1:TDRP.grd.nop
  % start and end of current patch in y/ydot vector and in cell
  % centre matrix
  ystart = TDRP.grd.ps(patchId);
  yend   = TDRP.grd.pe(patchId);
  ccstart = 1+(ystart-1)/TDRP.tdr.size;
  ccend   = yend/TDRP.tdr.size;
  
  % compute reaction terms in ydot(ystart:yend)
  if TDRP.tdr.haveReactionTerms
    % we have to compute some reaction term
    if (TDRP.tdr.FReac.vec)
      % we call the vectorized reaction computation
      ttmp=cputime;
      spaceIndices = [ystart/TDRP.tdr.size:(-1+(yend+1)/TDRP.tdr.size)];
      ydotreac = ProbFReac(t, ...
			   TDRP.grd.cellCentreMatrix(:,ccstart:ccend), ...
			   reshape(y(ystart:yend), TDRP.tdr.size, ...
				   ccend-ccstart+1), ...
			   TDRP.params);
      ydot(ystart:yend) = ydotreac(:);
      TIMER.tdrReacV = TIMER.tdrReacV + (cputime - ttmp);
      TIMER.tdrReacVC = TIMER.tdrReacVC + 1;
    else
      % we call the non-vectorized reaction computation
      ttmp=cputime;
      p = ystart-1;
      cc = ccstart-1;
      for cellId = 1:(ccend-ccstart+1)
	ydot(p+[1:TDRP.tdr.size]) = ...
	    ProbFReac(t, ...
		      TDRP.grd.cellCentreMatrix(:,cc+cellId), ...
		      y(p+[1:TDRP.tdr.size]), ...
		      TDRP.params);
	p = p + TDRP.tdr.size;
      end
      TIMER.tdrReac = TIMER.tdrReac + (cputime - ttmp);
      TIMER.tdrReacC = TIMER.tdrReacC + 1;
    end
  end % of if TDRP.tdr.haveReactionTerms

  % check if there is some transport or non-local term
  % if yes prepare for and do the transport discretisation
  if (TDRP.tdr.haveDiffusionTerms || ...
      TDRP.tdr.haveTaxisTerms || ...
      TDRP.tdr.haveNonLocalTerms) 
    % there are transport or non-local terms 
    % --> create TDRP.data{patchId} data structure as preparation for
    % transport/non-local discretisation 
    ttmp=cputime;
    tdrPrepTrans(t, y, patchId);
    TIMER.tdrPrepTrans = TIMER.tdrPrepTrans + (cputime - ttmp);
    TIMER.tdrPrepTransC = TIMER.tdrPrepTransC + 1;
  
    % discretise transport/non-local term on patch
    TDRP.data{patchId}.ComputedFaceData = false;

    % compute diffusion on current patch
    if TDRP.tdr.haveDiffusionTerms
      ttmp=cputime;
      tdrDiff(t, patchId);   % adds discretisation terms to
			     % TDRP.data.ydot{patchId} 
      TIMER.tdrDiff = TIMER.tdrDiff + (cputime - ttmp);
      TIMER.tdrDiffC = TIMER.tdrDiffC + 1; 
    end

    % compute taxis and/or non-local on current patch
    if (TDRP.tdr.haveTaxisTerms || TDRP.tdr.haveNonLocalTerms)
      ttmp=cputime;
      tdrTaxis(t, patchId); % adds discretisation terms to
			    % TDRP.data.ydot{patchId}  
      TIMER.tdrTaxis = TIMER.tdrTaxis + (cputime - ttmp);
      TIMER.tdrTaxisC = TIMER.tdrTaxisC + 1; 
    end
  
    % finalise transport discretisation by adding to ydot
    ttmp=cputime;
    if TDRP.grd.is1D
      ydot = ydot + reshape(TDRP.data{patchId}.ydot',length(ydot),1);
    else
      for s = 1:TDRP.tdr.size
	ydot(ystart+s-1:TDRP.tdr.size:yend) = ...
	    ydot(ystart+s-1:TDRP.tdr.size:yend) ...
	    + reshape(TDRP.data{patchId}.ydot(:,:,s), ...
		      TDRP.grd.nx(patchId)*TDRP.grd.ny(patchId), 1);
      end
    end
    TIMER.tdrFinTrans = TIMER.tdrFinTrans + (cputime - ttmp);
    TIMER.tdrFinTransC = TIMER.tdrFinTransC + 1;
  end % of 'if there is transport'
  
  % clear entry TDRP.data(patchId)
  TDRP.data{patchId} = [];
end % of for patchId = 1:TDRP.grd.nop

%%% HACK %%%
if norm(ydot,inf)==0 
  ydot(1)=ydot(1)+1e-14*(rand(1)-0.5);
end


return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tdrPrepTrans(t, y, patchId)
%
% function tdrPrepTrans(t, y, patchId)
%
% This function prepares the transport discretisation for patch
% patchId. It creates a struct array entry 'data(patchId)' of
% the global varaible TDRP and fills it with appropriate values as
% described below (bW = TDRP.grd.boundaryWidth): 
%   data(patchId).t    ... time t for which the data in
%                          data(patchId) is valid
%   data(patchId).ydot ... is a 3D matrix of 
% 			   dimensions grd.ny(patchId) times grd.nx(patchId) 
% 			   times tdr.size used to accumulate the 
% 			   contributions to ydot for each patch; 
% 			   initialised with zeros.
%   data(patchId).y    ... is a 3D matrix of dimensions 
% 			   grd.ny(patchId)+2*bW times grd.nx(patchId)+2*bW 
% 			   times tdr.size, initialised with solution 
% 			   values as provided in vector y and boundary 
% 			   values as derived from the boundary conditions
% 			   or values taken from neighbouring patches.
% 
% 
% This function calls external functions: 
%        ProbBCs()
% 
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : in tdr/tdrFdgl.m
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definition of Standard Codes               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defined boundaries
left   = 1;
right  = 2;
bottom = 3;
top    = 4;
% defined boundary condition types
None      = 0;
ZeroFlux  = 1;
Dirichlet = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End of Definition of Standard Codes        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear entry TDRP.data(patchId)
TDRP.data{patchId} = [];
% set current time
TDRP.data{patchId}.t = t;

% Get boundary width to use and determine ystart and yend such that
% y(ystart:yend) corresponds to the current patchId
bW     = TDRP.grd.boundaryWidth;
ystart = TDRP.grd.ps(patchId);
yend   = TDRP.grd.pe(patchId);

if (TDRP.grd.is1D)
  if (patchId ~= 1)
    error('... that should not have happened...');
  end
  TDRP.data{1}.y = NaN(TDRP.grd.ny(1) + 2*bW, ...
			     TDRP.tdr.size);
  TDRP.data{1}.y(bW+[1:TDRP.grd.ny(1)],:) = ...
      reshape(y,TDRP.tdr.size,TDRP.grd.ny(1))';
  TDRP.data{1}.ydot = zeros(TDRP.grd.ny(1), ...
				  TDRP.tdr.size);
  % fill the extended boundary in TDRP.data{1}.y
  %
  % left boundary 1D
  %
  ngbPatchId = TDRP.grd.ngb(1,left);
  if (ngbPatchId == 1)
  % periodic BCs
  TDRP.data{1}.y([1:bW],:) = ...
      TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[-bW+1:0],:);
      
  elseif (ngbPatchId ~= 0)
    error('That should not happen in 1D');
  else
    % compute boundary data
    yvals = TDRP.grd.y0(1);
    [BCtype,BCval] = ProbBCs(1, left, t, yvals, TDRP.params);
    % returns a length(yvals) times number of PDEs matrix BCType, where the
    % (i,j)th entry encodes the type of Boundary condition in cell i of
    % equation j. If BCtype(i,j) == Dirichlet then BCval(i,j) contains the
    % corresponding Dirichlet value, otherwise BCval(i,j) is unused.
    for eqnId=1:TDRP.tdr.size
      switch BCtype(1, eqnId)
       case None
	TDRP.data{1}.y([1:bW],eqnId) = NaN(bW,1);
       case Dirichlet
	TDRP.data{1}.y([1:bW],eqnId) = BCval(1, eqnId)*ones(bW,1);
       case ZeroFlux
	TDRP.data{1}.y([1:bW],eqnId) = ...
	    TDRP.data{1}.y(bW+[bW:-1:1],eqnId);
      otherwise
	error('tdrPrepTrans::left1D: unknown boundary condition type');
      end
    end
  end
  %
  % right boundary 1D
  %
  ngbPatchId = TDRP.grd.ngb(1,right);
  if (ngbPatchId == 1)
  % periodic BCs
  TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[1:bW],:) = ...
      TDRP.data{1}.y(bW+[1:bW],:);
  elseif (ngbPatchId ~= 0)
    error('That should not happen in 1D');
  else
    % compute boundary data
    yvals = TDRP.grd.y0(1)+TDRP.grd.dy(1)*TDRP.grd.ny(1);
    [BCtype,BCval] = ProbBCs(1, right, t, yvals, TDRP.params);
    % returns a length(yvals) times number of PDEs matrix BCType, where the
    % (i,j)th entry encodes the type of Boundary condition in cell i of
    % equation j. If BCtype(i,j) == Dirichlet then BCval(i,j) contains the
    % corresponding Dirichlet value, otherwise BCval(i,j) is unused.
    for eqnId=1:TDRP.tdr.size
      switch BCtype(1, eqnId)
       case None
	TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[1:bW],eqnId) = NaN(bW,1);
       case Dirichlet
	TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[1:bW],eqnId) = ...
	    BCval(1, eqnId)*ones(bW,1);
       case ZeroFlux
	TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[1:bW],eqnId) = ...
	    TDRP.data{1}.y(bW+TDRP.grd.ny(1)+[0:-1:-bW+1],eqnId);
      otherwise
	error('tdrPrepTrans::right1D: unknown boundary condition type');
      end
    end
  end
  % 1D stuff is finished here!
  return;
end

% Form extended 3D-matrix TDRP.data.y{patchId} for all concentrations 
% in current patch and initialise with NaN
TDRP.data{patchId}.y = NaN(TDRP.grd.ny(patchId) + 2*bW, ...
			   TDRP.grd.nx(patchId) + 2*bW, ...
			   TDRP.tdr.size);
% Fill the centre part with the given y-values
for s = 1:TDRP.tdr.size
  TDRP.data{patchId}.y(bW + [1:TDRP.grd.ny(patchId)], ...
		       bW + [1:TDRP.grd.nx(patchId)], ...
		       s) ...
      = reshape(y(ystart+s-1:TDRP.tdr.size:yend), ...
		TDRP.grd.ny(patchId), TDRP.grd.nx(patchId));
end

% Create 3D-matrix TDRP.data{patchId}.ydot to accumulate the
% transport discretisation and initialise with zeros
TDRP.data{patchId}.ydot = zeros(TDRP.grd.ny(patchId), ...
				TDRP.grd.nx(patchId), ...
				TDRP.tdr.size);

% fill the extended boundary in patch patchId 
%
% left boundary
%
ngbPatchId = TDRP.grd.ngb(patchId,left);
if ( ngbPatchId ~= 0)
  % copy from other patch
  rows = [1:TDRP.grd.ny(patchId)];
  cols = TDRP.grd.nx(ngbPatchId)+[-bW+1:0];
  spec = [1:TDRP.tdr.size];
  TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)],[1:bW],:) = ...
      tdrGetSolutionValues(y, ngbPatchId, rows, cols, spec);
else
  % compute boundary data
  xvals = TDRP.grd.x0(patchId) * ones(TDRP.grd.ny(patchId), 1);
  yvals = TDRP.grd.yc{patchId};
  [BCtype,BCval] = ProbBCs(patchId, left, t, xvals, yvals, TDRP.params);
  % returns a length(xvals) * number of PDEs matrix BCType, where the
  % (i,j)th entry encodes the type of Boundary condition in cell i of
  % equation j. If BCtype(i,j) == Dirichlet then BCval(i,j) contains the
  % corresponding Dircihlet value, otherwise BCval(i,j) is unused.
  for eqnId=1:TDRP.tdr.size
    BCNoneIndex      = find(BCtype(:, eqnId) == None);
    BCDirichletIndex = find(BCtype(:, eqnId) == Dirichlet);
    BCZeroFluxIndex  = find(BCtype(:, eqnId) == ZeroFlux);
    if (length(BCDirichletIndex)+length(BCZeroFluxIndex) ...
	+ length(BCNoneIndex) ~= TDRP.grd.ny(patchId))
      error('tdrPrepTrans::left: unknown boundary condition type');
    end
    % deal with ZeroFlux
    if (length(BCZeroFluxIndex))
      TDRP.data{patchId}.y      (bW+BCZeroFluxIndex,    [1:bW]   , eqnId) ...
	  = TDRP.data{patchId}.y(bW+BCZeroFluxIndex, bW+[bW:-1:1], eqnId);
    end
    % deal with Dirichlet
    if (length(BCDirichletIndex))
      TDRP.data{patchId}.y      (bW+BCDirichletIndex,   [1:bW]   , eqnId) ...
	  = BCval(BCDirichletIndex, eqnId)*ones(1,bW);
    end
    % deal with None
    if (length(BCNoneIndex))
      TDRP.data{patchId}.y      (bW+BCNoneIndex,        [1:bW]   , eqnId) ...
	  = NaN * ones(length(BCNoneIndex), bW);
    end
  end
end
  
%
% right boundary
%
ngbPatchId = TDRP.grd.ngb(patchId,right);
if (ngbPatchId ~= 0)
  % copy from other patch
  rows = [1:TDRP.grd.ny(patchId)];
  cols = [1:bW];
  spec = [1:TDRP.tdr.size];
  TDRP.data{patchId}.y(bW+[1:TDRP.grd.ny(patchId)],end+[-bW+1:0],:) = ...
      tdrGetSolutionValues(y, ngbPatchId, rows, cols, spec);
else
  % compute boundary data
  xvals = (TDRP.grd.x0(patchId) + TDRP.grd.nx(patchId)*TDRP.grd.dx(patchId))...
	  * ones(TDRP.grd.ny(patchId), 1);
  yvals = TDRP.grd.yc{patchId};
  [BCtype,BCval] = ProbBCs(patchId, right, t, xvals, yvals, TDRP.params);
  for eqnId=1:TDRP.tdr.size
    BCNoneIndex      = find(BCtype(:, eqnId) == None);
    BCDirichletIndex = find(BCtype(:, eqnId) == Dirichlet);
    BCZeroFluxIndex  = find(BCtype(:, eqnId) == ZeroFlux);
    if (length(BCDirichletIndex)+length(BCZeroFluxIndex) ...
	+ length(BCNoneIndex) ~= TDRP.grd.ny(patchId))
      error('tdrPrepTrans::right: unknown boundary condition type');
    end
    % deal with ZeroFlux
    if (length(BCZeroFluxIndex))
      TDRP.data{patchId}.y      (bW+BCZeroFluxIndex,    end+[-bW+1:0]   , eqnId) ...
	  = TDRP.data{patchId}.y(bW+BCZeroFluxIndex, end-bW+[0:-1:-bW+1], eqnId);
    end
    % deal with Dirichlet
    if (length(BCDirichletIndex))
      TDRP.data{patchId}.y      (bW+BCDirichletIndex,   end+[-bW+1:0]   , eqnId) ...
	  = BCval(BCDirichletIndex, eqnId)*ones(1,bW);
    end
    % deal with None
    if (length(BCNoneIndex))
      TDRP.data{patchId}.y      (bW+BCNoneIndex,        end+[-bW+1:0]   , eqnId) ...
	  = NaN * ones(length(BCNoneIndex), bW);
    end
  end
end
  
%
% bottom boundary
%
ngbPatchId = TDRP.grd.ngb(patchId, bottom);
if (ngbPatchId ~= 0)
  % copy from other patch
  rows = TDRP.grd.ny(ngbPatchId)+[-bW+1:0];
  cols = [1:TDRP.grd.nx(patchId)];
  spec = [1:TDRP.tdr.size];
  TDRP.data{patchId}.y([1:bW], bW+[1:TDRP.grd.nx(patchId)], :) = ...
      tdrGetSolutionValues(y, ngbPatchId, rows, cols, spec);
else
  % compute boundary data
  xvals = TDRP.grd.xc{patchId};
  yvals = TDRP.grd.y0(patchId) * ones(TDRP.grd.nx(patchId), 1);
  [BCtype,BCval] = ProbBCs(patchId, bottom, t, xvals, yvals, TDRP.params);
  for eqnId=1:TDRP.tdr.size
    BCNoneIndex      = find(BCtype(:, eqnId) == None);
    BCDirichletIndex = find(BCtype(:, eqnId) == Dirichlet);
    BCZeroFluxIndex  = find(BCtype(:, eqnId) == ZeroFlux);
    if (length(BCDirichletIndex)+length(BCZeroFluxIndex) ...
	+ length(BCNoneIndex) ~= TDRP.grd.nx(patchId))
      error('tdrPrepTrans::bottom: unknown boundary condition type');
    end
    % deal with ZeroFlux
    if (length(BCZeroFluxIndex))
      TDRP.data{patchId}.y      (   [1:bW]   , bW+BCZeroFluxIndex, eqnId) ...
	  = TDRP.data{patchId}.y(bW+[bW:-1:1], bW+BCZeroFluxIndex, eqnId);
    end
    % deal with Dirichlet
    if (length(BCDirichletIndex))
      TDRP.data{patchId}.y      (   [1:bW]   ,bW+BCDirichletIndex, eqnId) ...
	  = ones(bW, 1) * BCval(BCDirichletIndex, eqnId)';
    end
    % deal with None
    if (length(BCNoneIndex))
      TDRP.data{patchId}.y      (   [1:bW]   , bW+BCNoneIndex    , eqnId) ...
	  = NaN * ones(bW, length(BCNoneIndex));
    end
  end
end
  
%
% top boundary
%
ngbPatchId = TDRP.grd.ngb(patchId,top);
if (ngbPatchId ~= 0)
  % copy from other patch
  rows = [1:bW];
  cols = [1:TDRP.grd.nx(patchId)];
  spec = [1:TDRP.tdr.size];
  TDRP.data{patchId}.y(end+[-bW+1:0], bW+[1:TDRP.grd.nx(patchId)], :) = ...
      tdrGetSolutionValues(y, ngbPatchId, rows, cols, spec);
else
  % compute boundary data
  xvals = TDRP.grd.xc{patchId};
  yvals = (TDRP.grd.y0(patchId) + TDRP.grd.dy(patchId) * TDRP.grd.ny(patchId))* ...
	  ones(TDRP.grd.nx(patchId), 1);
  [BCtype,BCval] = ProbBCs(patchId, top, t, xvals, yvals, TDRP.params);
  for eqnId=1:TDRP.tdr.size
    BCNoneIndex      = find(BCtype(:, eqnId) == None);
    BCDirichletIndex = find(BCtype(:, eqnId) == Dirichlet);
    BCZeroFluxIndex  = find(BCtype(:, eqnId) == ZeroFlux);
    if (length(BCDirichletIndex)+length(BCZeroFluxIndex) ...
      + length(BCNoneIndex) ~= TDRP.grd.nx(patchId))
      error('tdrPrepTrans::top: unknown boundary condition type');
    end
    % deal with ZeroFlux
    if (length(BCZeroFluxIndex))
      TDRP.data{patchId}.y      (end+[-bW+1:0]      , bW+BCZeroFluxIndex, eqnId) ...
	  = TDRP.data{patchId}.y(end-bW+[0:-1:-bW+1], bW+BCZeroFluxIndex, eqnId);
    end
    % deal with Dirichlet
    if (length(BCDirichletIndex))
      TDRP.data{patchId}.y      (end+[-bW+1:0]      ,bW+BCDirichletIndex, eqnId) ...
	  = ones(bW, 1) * BCval(BCDirichletIndex, eqnId)';
    end
    % deal with None
    if (length(BCNoneIndex))
      TDRP.data{patchId}.y      (end+[-bW+1:0]      , bW+BCNoneIndex    , eqnId) ...
	  = NaN * ones(bW, length(BCNoneIndex));
    end
  end
end

% 2D stuff is finished now  
return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
