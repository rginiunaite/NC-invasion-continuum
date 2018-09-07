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
% The model consists of 1 PDEs implemented in the following order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Indices of PDEs                %%%%
 params.eq.c1  = 1;
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
params.domainlength = setValue('domainlength', 1, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.domainlength = ' ...
      num2str(params.domainlength) '.']);
params.gridCells = setValue('gridCells', 400, inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.gridCells = ' ...
      num2str(params.gridCells) ' grid cells per unit length.']);
params.BCs = setValue('BCs', 'vz', inputStruct, 'params');
disp(['ProbGetParams()::Info: Selected params.BCs = ' ...
      num2str(params.BCs) '.']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scaling parameters                      %%%
%params.scal.T    = 10000.0;   % [seconds]
%%% Scaling parameters                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 1. initialise struct tdr
%
tdr.versionString = '1.6pre-AdhPack-1.2';
tdr.size   = 1; % number of species (PDEs)
tdr.tvec   = [0:0.01:1 2:2:10]; % initial time and output times.


%%% CODES %%%%
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
tdr.FReac.vec = true;  % true if ProbFReac.m can be called vectorized

% describe transport function ProbFTrans.m
% if only Zero entries then ProbFTrans.m is never called
tdr.FTrans.depends = Zero * ones(tdr.size, tdr.size);
tdr.FTrans.depends(params.eq.c1, params.eq.c1) = Const; % constant diffusion of c1
tdr.FTrans.vec = true; % true if ProbFTrans.m can be called vectorized

% describe nonlocal function ProbFNonLocal.m
tdr.FNonLocal.depends = Zero * ones(tdr.size, 1);
tdr.FNonLocal.depends(params.eq.c1) = DependsU;
tdr.FNonLocal.R          = 0.1; 
tdr.FNonLocal.OmegaStr = [...
    '@(r)(1/(2*' num2str(tdr.FNonLocal.R, '%25.16d') ') * ' ...
    '(r <= ' num2str(tdr.FNonLocal.R,'%25.16d') ') )'];
tdr.FNonLocal.BCs = params.BCs;
tdr.FNonLocal

%
% 2. initialise struct params
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true
  %%% Parameters for no aggregation (higher initial cell density)
  params.IC      = 0.5;  % initial cell density
  params.c1_D    = 0.01;  % cell random motility
  params.c_Sc1c1 = 3;   % cell-cell adhesion coefficient 
else
  %%% Parameters for aggregation (lower initial cell density)
  params.IC      = 0.1;  % initial cell density
  params.c1_D    = 0.01; % cell random motility
  params.c_Sc1c1 = 10;   % cell-cell adhesion coefficient 
end
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
grd.is1D = false;
grd.x0  = 0.0;
grd.y0  = 0.0;
xfac=1;  % factor for different length in x-direction
grd.nx  = round(xfac*params.domainlength*params.gridCells);
grd.ny  = round(1*params.domainlength*params.gridCells);
grd.dx  = (xfac*params.domainlength) / grd.nx;
grd.dy  = (1*params.domainlength) / grd.ny;
switch params.BCs
 case {'pppp', 'pp'}
  grd.ngb = [ 1 1 1 1];
 case {'zzzz'}
  grd.ngb = [ 0 0 0 0];
 otherwise
  error('Unsupported BCs');
end
grd
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
