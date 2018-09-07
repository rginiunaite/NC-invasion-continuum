function [y0, tspan] = tdrInit(varargin)
%
% function [y0, tspan] = tdrInit(varargin);
%
% This function initialises a TDR-system for subsequent time integration
% by setting up a global struct TDRP with the data provided by a call to
% problem function ProbGetParams(). All input arguments are passed 
% unchanged to ProbGetParams(). Some extra fields are added to substruct
% TDRP.grd . The initial data is computed by calls to problem function 
% ProbFy0() and returned in standard order in vector y0. Also returned is
% the vector of output times tspan (returned as tdr.tvec by 
% ProbGetParams()). Both, y0 and tspan, are required for a subsequent 
% call to a time integration subroutine).
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdr/tdrInit.m
%* Date created  : 2005, July 13
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.1 (AG, 2008 09 08)
%*                     - Introduced TDRP.tdr.have{Reaction,Diffusion, 
%*                       Taxis,NonLocal}Terms in the TDRP data structure 
%*                       (defined in this function)
%*                     - Define TDRP.grd.boundaryWidth here (was in
%*                       TDRP.data before and defined elsewhere).
%*                     - added check for grid suitability in case
%*                       of the presence of non-local terms. If all
%*                       ok then setup required data for their
%*                       evaluation in struct TDRP.tdr.FNonLocal
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


% clear TDRP if already present
clear global TDRP

% now create it again as global
global TDRP;           % global struct for TDR Problem data


% First get the basic grid data, TDR data, and problem specific data 
% by a call to the problem function ProbGetParams()
[TDRP.grd, TDRP.tdr, TDRP.params] = ProbGetParams(varargin{:});

% Check whether model is for this version of the mTDR system
if ~strcmp(TDRP.tdr.versionString, tdrVersion)
  errmsg = ['tdrInit::model is for mTDR version ' TDRP.tdr.versionString ...
	    ' but current implementation is ' tdrVersion];
  error(errmsg);
else
  disp(['tdrInit::running current model with mTDR version ' ...
	tdrVersion]);
end

% Check that only allowed subfields are present in the user data
% structures 
%
% TODO

%
% Set optional grd parameters to default values if not present
%
if ~isfield(TDRP.grd, 'is1D')
  TDRP.grd.is1D = false;
end
if ~isfield(TDRP.grd, 'isAxisymmetric')
  TDRP.grd.isAxisymmetric = false;
end

%
% check what terms we have in the PDEs
%
TDRP.tdr.haveReactionTerms = false;
if sum(TDRP.tdr.FReac.depends)
  TDRP.tdr.haveReactionTerms = true;
end

TDRP.tdr.haveDiffusionTerms = false;
if sum(diag(TDRP.tdr.FTrans.depends))
  TDRP.tdr.haveDiffusionTerms = true;
end

TDRP.tdr.haveTaxisTerms = false;
if (sum(TDRP.tdr.FTrans.depends(:)) ~= ...
    sum(diag(TDRP.tdr.FTrans.depends)))
  TDRP.tdr.haveTaxisTerms = true;
end

TDRP.tdr.haveNonLocalTerms = false;
if isfield(TDRP.tdr, 'FNonLocal')
  if isfield(TDRP.tdr.FNonLocal, 'depends')
    if sum(TDRP.tdr.FNonLocal.depends(:))
      TDRP.tdr.haveNonLocalTerms = true;
    end
  end
end


%
% Set boundary width for patches
%
TDRP.grd.boundaryWidth = 0;
% check if there is diffusion
if TDRP.tdr.haveDiffusionTerms
  TDRP.grd.boundaryWidth = max(1,TDRP.grd.boundaryWidth);
end
% check if there is taxis
if TDRP.tdr.haveTaxisTerms
  TDRP.grd.boundaryWidth = max(2,TDRP.grd.boundaryWidth);
end
% check if there are non-local terms
if TDRP.tdr.haveNonLocalTerms
  TDRP.grd.boundaryWidth = max(2,TDRP.grd.boundaryWidth);
end


%
% Setup additional grd data (1D/2D) and check if we have a problem
% in 1D space 
%
if ~TDRP.grd.is1D % we have a 2D problem
  % add vectors ps (PatchStart)  and pe (PatchEnd) to TDRP.grd: 
  % TDRP.grd.ps(patchId) is the index in the solution vector of the 
  % first solution component in patch patchId
  % TDRP.grd.pe(patchId) is the index in the solution vector of the 
  % last solution component in patch patchId
  gridsize = 0;  % total number of grid cells in domain
  molsize = 0;   % dimension of the MOL-ODE system
  TDRP.grd.ps(1) = 1;
  for patchId = 1:TDRP.grd.nop
    gridsize = gridsize + TDRP.grd.nx(patchId) * TDRP.grd.ny(patchId);
    molsize = gridsize * TDRP.tdr.size;
    TDRP.grd.pe(patchId) = molsize;
    if (patchId < TDRP.grd.nop)
      TDRP.grd.ps(patchId+1) = molsize + 1;
    end
  end
  TDRP.grd.gridsize = gridsize;
  % Create vectors with x- and y-cell centres for each patch and
  % add to TDRP.grd 
  for patchId = 1:TDRP.grd.nop
    TDRP.grd.xc{patchId} = TDRP.grd.x0(patchId) ...
	+ TDRP.grd.dx(patchId)*([1:TDRP.grd.nx(patchId)]-0.5);
    TDRP.grd.yc{patchId} = TDRP.grd.y0(patchId) ...
	+ TDRP.grd.dy(patchId)*([1:TDRP.grd.ny(patchId)]-0.5);
  end
  % Create matrix with cell centres (x-coordinates in first row,
  % y-coordinates in second row) and add to TDRP.grd. Within each 
  % patch we store y-value changing faster than x-value.
  TDRP.grd.cellCentreMatrix = zeros(2, gridsize);
  p = 1; % pointer in TDRP.grd.cellCentreMatrix
  for patchId = 1:TDRP.grd.nop
    for j = 1:TDRP.grd.nx(patchId)
      for i = 1:TDRP.grd.ny(patchId)
	TDRP.grd.cellCentreMatrix(:, p) = ...
	    [TDRP.grd.xc{patchId}(j); TDRP.grd.yc{patchId}(i)];
	p = p + 1;
      end
    end
  end    
else % we have a 1D problem
  if TDRP.grd.isAxisymmetric
    error(['tdrInit:: Axisymmetry is not (yet) supported for 1D' ...
	   ' domains.']); 
  end
  if (TDRP.grd.nop ~= 1)
    error(['tdrInit:: 1D domains allow for one patch only.']);
  end
  if (size(TDRP.grd.ngb) ~= [1 2])
    error(['tdrInit:: 1D domains require 1 x 2 neighbour vector.']);
  end
  % setup additional grd data
  TDRP.grd.ps(1) = 1;
  TDRP.grd.pe(1) = TDRP.tdr.size*TDRP.grd.ny(1);
  TDRP.grd.gridsize = TDRP.grd.ny(1);
  TDRP.grd.yc{1} = TDRP.grd.y0(1) + TDRP.grd.dy(1)*([1:TDRP.grd.ny(1)]-0.5);
  TDRP.grd.cellCentreMatrix = TDRP.grd.yc{1};
end

%
% Check grid data for suitability with non-local term and setup
% non-local term required data
%
if TDRP.tdr.haveNonLocalTerms
  % we have non-local terms, so make sure the grid is suitable for 
  % them 
  if TDRP.grd.isAxisymmetric
    error(['tdrInit:: The presence of non-local terms does not allow' ...
	   ' for axisymmetric domains']);
  end
  if (TDRP.grd.nop ~= 1)
    error(['tdrInit:: The presence of non-local terms allows only' ...
	   ' for one patch.']);
  end
  if TDRP.grd.is1D
    if (TDRP.grd.ngb ~= [1 1]) & (TDRP.grd.ngb ~= [0 0])
      error(['tdrInit:: The presence of non-local terms allows only' ...
	     ' for periodic boundary conditions or BCs at both end.']);
    end
    % setup of masks for non-local term evaluation
    % ... here we assume that the non-local term in each equation has
    % the same form except for the function g
    disp(['tdrInit: Starting setup of mask for non-local term' ...
	  ' evaluation...']); 
    if (size(TDRP.tdr.FNonLocal.R) == [1,1])
      disp('tdrInit: only one type of nonlocal term.');
      TDRP.tdr.FNonLocal.h = TDRP.grd.dy(1);
      TDRP.tdr.FNonLocal.N2 = TDRP.grd.ny(1);
      TDRP.tdr.FNonLocal.RuleId = 6;
      TDRP.tdr.FNonLocal.NR = max(1000,...
				  round(TDRP.tdr.FNonLocal.N2 ...
					*TDRP.tdr.FNonLocal.R)...
				  );
      mask = ...
	  setupIntegralRule1D_weights(TDRP.tdr.FNonLocal.R, ...
				      TDRP.tdr.FNonLocal.h, ...
				      eval(TDRP.tdr.FNonLocal.OmegaStr), ...
				      TDRP.tdr.FNonLocal.RuleId, ...
				      TDRP.tdr.FNonLocal.NR);
      TDRP.tdr.FNonLocal.mask = ...
	  setupIntegralRule1D_BCs(mask, TDRP.tdr.FNonLocal.N2, ...
					TDRP.tdr.FNonLocal.BCs);
    else
      disp('tdrInit: multiple types of nonlocal terms.');
      TDRP.tdr.FNonLocal.h = TDRP.grd.dy(1);
      TDRP.tdr.FNonLocal.N2 = TDRP.grd.ny(1);
      TDRP.tdr.FNonLocal.RuleId = 6;
      Rmax = max(TDRP.tdr.FNonLocal.R(:));
      TDRP.tdr.FNonLocal.NR = max(1000,...
				  round(TDRP.tdr.FNonLocal.N2*Rmax));
      if ~strcmp(TDRP.tdr.FNonLocal.BCs,'pp')
	error(['tdrInit: multiple nonlocal terms only with periodic' ...
	       ' BCs.'])
      end
      TDRP.tdr.FNonLocal.mask = cell(size(TDRP.tdr.FNonLocal.R));
      for ii=1:size(TDRP.tdr.FNonLocal.R,1)
	for jj=1:size(TDRP.tdr.FNonLocal.R,2)
	  if (TDRP.tdr.FNonLocal.R(ii,jj)>0)
	    mask = ...
		setupIntegralRule1D_weights(TDRP.tdr.FNonLocal.R(ii,jj), ...
			TDRP.tdr.FNonLocal.h, ...
			eval(TDRP.tdr.FNonLocal.OmegaStr{ii,jj}), ...
			TDRP.tdr.FNonLocal.RuleId, ...
			TDRP.tdr.FNonLocal.NR);
	    TDRP.tdr.FNonLocal.mask{ii,jj} = ...
		setupIntegralRule1D_BCs(mask, TDRP.tdr.FNonLocal.N2, ...
					      TDRP.tdr.FNonLocal ...
					      .BCs);
	  else
	    TDRP.tdr.FNonLocal.mask{ii,jj} = NaN;
	  end
	end
      end
    end
    disp(['tdrInit: Completed setup of mask for non-local term' ...
	  ' evaluation.']); 
  else % the grid is 2D
    if (TDRP.grd.ngb ~= [1 1 1 1]) & (TDRP.grd.ngb ~= [0 0 0 0])
      error(['tdrInit:: The presence of non-local terms allows only' ...
	     ' for periodic boundary conditions or one patch with' ...
	     ' BCs on all sides.']);
    end
    if (abs(TDRP.grd.dx(1)-TDRP.grd.dy(1))>=10*eps)
      error(['tdrInit:: The presence of non-local terms allows only' ...
	     ' equal spatial grid width in both directions.']);
    else
      if (TDRP.grd.dx(1) ~= TDRP.grd.dy(1))
	disp(['tdrInit:: grd.dx(1) and grd.dy(1) differ by less than' ...
	      ' 10*eps. Setting them equal.'])
      end
    end
    if (size(TDRP.tdr.FNonLocal.R) ~= [1,1])
      error('tdrInit: multiple nonlocal terms not supported in 2D.');
    end
    % setup of masks for non-local term evaluation
    % ... here we assume that the non-local term in each equation has
    % the same form except for the function g
    disp(['tdrInit: Starting setup of masks for non-local term' ...
	  ' evaluation...']); 
    TDRP.tdr.FNonLocal.h = TDRP.grd.dx(1);
    TDRP.tdr.FNonLocal.N1 = TDRP.grd.nx(1);
    TDRP.tdr.FNonLocal.N2 = TDRP.grd.ny(1);
    TDRP.tdr.FNonLocal.RuleId = 1;
    TDRP.tdr.FNonLocal.NR = 100;
    TDRP.tdr.FNonLocal.weightTol = 5e-5;
    TDRP.tdr.FNonLocal.weightTol = 5e-4;
    
    switch TDRP.tdr.FNonLocal.BCs
      case 'pp'
       TDRP.tdr.FNonLocal.mask.x1Dir = ...
	   setupIntegralRule2D_old('x1', TDRP.tdr.FNonLocal.R, ...
					 TDRP.tdr.FNonLocal.h, ...
					 TDRP.tdr.FNonLocal.N1, ...
					 TDRP.tdr.FNonLocal.N2, ...
					 eval(TDRP.tdr.FNonLocal.OmegaStr), ...
					 TDRP.tdr.FNonLocal.RuleId, ...
					 TDRP.tdr.FNonLocal.NR, ...
					 TDRP.tdr.FNonLocal.weightTol);
       TDRP.tdr.FNonLocal.mask.x2Dir = ...
	   setupIntegralRule2D_old('x2', TDRP.tdr.FNonLocal.R, ...
					 TDRP.tdr.FNonLocal.h, ...
					 TDRP.tdr.FNonLocal.N1, ...
					 TDRP.tdr.FNonLocal.N2, ...
					 eval(TDRP.tdr.FNonLocal.OmegaStr), ...
					 TDRP.tdr.FNonLocal.RuleId, ...
					 TDRP.tdr.FNonLocal.mask.x1Dir.NR, ...
					 inf);
     case {'pppp', 'zzzz'}
      clear TDRP.tdr.FNonLocal.mask
      TDRP.tdr.FNonLocal.mask.BCs = TDRP.tdr.FNonLocal.BCs;
      nonlocal = setupIntegralRule2D_weights('x1', ...
					     TDRP.tdr.FNonLocal.R, ...
					     TDRP.tdr.FNonLocal.h, ...
					     eval(TDRP.tdr.FNonLocal.OmegaStr), ...
					     TDRP.tdr.FNonLocal.RuleId, ...
					     TDRP.tdr.FNonLocal.NR, ...
					     TDRP.tdr.FNonLocal.weightTol);
      TDRP.tdr.FNonLocal.mask.x1Dir = ...
	  setupIntegralRule2D_BCs(nonlocal, ...
				  TDRP.tdr.FNonLocal.N1, ...
				  TDRP.tdr.FNonLocal.N2, ...
				  TDRP.tdr.FNonLocal.BCs);
      nonlocal = setupIntegralRule2D_weights('x2', ...
					     TDRP.tdr.FNonLocal.R, ...
					     TDRP.tdr.FNonLocal.h, ...
					     eval(TDRP.tdr.FNonLocal.OmegaStr), ...
					     TDRP.tdr.FNonLocal.RuleId, ...
					     TDRP.tdr.FNonLocal.mask.x1Dir.NR, ...
					     inf);
      TDRP.tdr.FNonLocal.mask.x2Dir = ...
	  setupIntegralRule2D_BCs(nonlocal, ...
				  TDRP.tdr.FNonLocal.N1, ...
				  TDRP.tdr.FNonLocal.N2, ...
				  TDRP.tdr.FNonLocal.BCs);
     otherwise
      error(['tdrInit:: Unsupported value of tdr.FNonLocal.BCs']);
    end
    disp(['tdrInit: Completed setup of masks for non-local term' ...
	  ' evaluation.']); 
  end % if TDRP.grd.is1D
end % TDRP.tdr.haveNonLocalTerms

% y0 initialisieren
y0 = zeros(TDRP.grd.gridsize * TDRP.tdr.size,1);
p = 1;
for cellId = 1:TDRP.grd.gridsize
  y0(p:p + TDRP.tdr.size - 1) = ...
      ProbFy0(TDRP.grd.cellCentreMatrix(:, cellId), TDRP.params);
  p = p + TDRP.tdr.size;
end

% tspan initialisieren
tspan = TDRP.tdr.tvec;

% create TDRP.data array 
TDRP.data=cell(TDRP.tdr.size,1);

return
% end of function
