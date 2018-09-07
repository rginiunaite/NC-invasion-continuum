%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [grd, tdr, params] = ProbGetParams(varargin)
%
% Input arguments
%   varargin      ... 
% Output arguments
%   grd           ... information on grid composition and size
%   tdr           ... information on the tdr system
%   params        ... parameter values of the problem
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grd, tdr, params] = ProbGetParams(varargin)

grd    = struct();
tdr    = struct();
params = struct();

%
% The model consists of 2 PDEs implemented in the following order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indices of PDEs                %%%%
params.eq.u  = 1;
params.eq.v  = 2;
params.eq.c  = 3;
%%% End of Indices of PDEs         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% process input parameters
switch length(varargin)
  case 0
    inputStruct.filled = false;
  case 1
    inputStruct = varargin{1};
    inputStruct.filled = true;
  otherwise
    error('Too many arguments in call to ProbGetParams().');
end


%
% initialise struct params
%
params.gridCells = setValue('gridCells', 400, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.gridCells = ' ...
      num2str(params.gridCells) ' grid cells per unit length.']);

%
% 1. initialise struct tdr
%
tdr.versionString = '1.6pre-AdhPack-1.2';
tdr.size   = 3; % number of species (PDEs)
tdr.tvec   = [0:0.001:0.1]; % initial, output and final time.

%%% CODES for describing a function f %%%%
Zero     = 0;              % f \equiv 0
Const    = 1;              % f \equiv const
DependsT = 2;              % f = f(t)
DependsS = 4;              % f = f(x,y)
DependsU = 8;              % f = f(u)
DependsTS = DependsT + DependsS; %  f = f(t, x, y)
DependsTU = DependsT + DependsU; %  f = f(t, u)
DependsSU = DependsS + DependsU; %  f = f(x, y, u)
DependsTSU = DependsT + DependsS + DependsU; %  f = f(t, x, y, u)

% describe reaction function ProbFReac.m 
% if only Zero entries then ProbFReac.m is never called
tdr.FReac.depends = Zero * ones(tdr.size, 1); 
tdr.FReac.depends(params.eq.u) = [ DependsTU ];
tdr.FReac.depends(params.eq.v) = [ DependsTU ];
tdr.FReac.depends(params.eq.c) = [ DependsTU ];
tdr.FReac.vec = true;  % true if ProbFReac.m can be called vectorized

% describe transport function ProbFTrans.m
% if only Zero entries then ProbFTrans.m is never called
tdr.FTrans.depends = Zero * ones(tdr.size, tdr.size);
tdr.FTrans.depends(params.eq.u, params.eq.u) = DependsT; % diff. coeff. of u
tdr.FTrans.depends(params.eq.v, params.eq.v) = DependsT; % diff. coeff. of v
tdr.FTrans.depends(params.eq.c, params.eq.c) = DependsT; % diff. coeff. of c
tdr.FTrans.depends(params.eq.u, params.eq.c) = DependsTU; % chemo of leaders to VEGF
tdr.FTrans.depends(params.eq.v, params.eq.u) = DependsTU; % guidance of followers by leader gradient
tdr.FTrans.vec = true; % true if ProbFTrans.m can be called vectorized

% describe nonlocal function ProbFNonLocal.m
% if only Zero entries then ProbFNonLocal.m is never called
% (ProbFNonLocal.m encodes the function g under the integral)
tdr.FNonLocal.depends = Zero * ones(tdr.size, 1);

% describe nonlocal function ProbFNonLocal2.m
% if only Zero entries then ProbFNonLocal2.m is never called and
% value 1 is assumed.
% (ProbFNonLocal.m encodes the function p multiplying the integral)
tdr.FNonLocal.depends2 = Zero * ones(tdr.size, 1);
%tdr.FNonLocal.depends2(params.eq.n1) = DependsU;
%tdr.FNonLocal.depends2(params.eq.n2) = DependsU;

% describe sensing radii
%tdr.FNonLocal.R          = zeros(tdr.size); 
%tdr.FNonLocal.R(params.eq.n1,params.eq.n1) = 0.1; % for equ. 1
%tdr.FNonLocal.R(params.eq.n1,params.eq.n2) = 0.1; % for equ. 1
%tdr.FNonLocal.R(params.eq.n2,params.eq.n1) = 0.1; % for equ. 2
%tdr.FNonLocal.R(params.eq.n2,params.eq.n2) = 0.1; % for equ. 2

% describe Omega function (the same for all)
%tdr.FNonLocal.OmegaStr = cell(tdr.size);
%tdr.FNonLocal.OmegaStr{params.eq.n1,params.eq.n1} = ['@(r)( (r <= ' ...
%    num2str(tdr.FNonLocal.R(params.eq.n1,params.eq.n1),'%25.16d') ...
%		    ') )'];
%tdr.FNonLocal.OmegaStr{params.eq.n1,params.eq.n2} = ['@(r)( (r <= ' ...
%    num2str(tdr.FNonLocal.R(params.eq.n1,params.eq.n2),'%25.16d') ...
%		    ') )'];
%tdr.FNonLocal.OmegaStr{params.eq.n2,params.eq.n1} = ['@(r)( (r <= ' ...
%    num2str(tdr.FNonLocal.R(params.eq.n2,params.eq.n1),'%25.16d') ...
%		    ') )'];
%tdr.FNonLocal.OmegaStr{params.eq.n2,params.eq.n2} = ['@(r)( (r <= ' ...
%    num2str(tdr.FNonLocal.R(params.eq.n2,params.eq.n2),'%25.16d') ...
%		    ') )'];

% specify boundary condition (only periodic supported)
%tdr.FNonLocal.BCs = params.BCs;
%tdr.FNonLocal

%
% 2. initialise struct params
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters of the initial conditions    %%%
params.IC.c       = 1;   % initial chemoattractant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters of the differential equation %%%
% reaction terms
params.reac_production_kappa = 1e-3;
params.reac_degradation_alpha = 10.0;
% cell random motility
params.u_D = 0.0001
params.v_D = 0.01;
% VEGF diffusion
params.c_D = 0.1;
% leader chemotaxis
params.selectChemotaxisLeader = 'receptor';
%params.selectChemotaxisLeader = 'simple';
params.chemotaxisLeader_k = 50;
params.chemotaxisLeader_beta = 1;
% follower chemotaxis
params.selectChemotaxisFollower = 'receptor';
%params.selectChemotaxisFollower = 'simple';
params.chemotaxisFollower_k = 50;
params.chemotaxisFollower_beta = 1;
% params domain growth
params.Linfty = 86.76;
params.a = 0.2888;
params.ts = 12.77;
params.L0 = 30;
params.k0 = 29.12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% params struct completed
%

%
% 3. define the spatial domain via its patches and store data in grd struct
%

% We use a 1-patch
grd.nop = 1;             % Number Of Patches
grd.isAxiSymmetric = false;
grd.is1D = true;
grd.y0  = 0.0;
grd.ny  = round(1*params.gridCells);
grd.dy  = 1/grd.ny;
grd.ngb = [ 0 0 ];
%
% grd struct completed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = setValue(fieldName, defaultValue, inputStruct, subStructName)
  val = defaultValue;
  if isfield(inputStruct, subStructName)
    tmpStruct = getfield(inputStruct, subStructName);
    if isfield(tmpStruct, fieldName)
      val = getfield(tmpStruct, fieldName);
    end
  end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
