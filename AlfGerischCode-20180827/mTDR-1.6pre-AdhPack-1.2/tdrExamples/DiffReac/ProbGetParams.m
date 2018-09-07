function [grd, tdr, params] = ProbGetParams(varargin)
% function [grd, tdr, params] = ProbGetParams(varargin)
%
% This function is part of the DiffReac TDR model.
%
% The model consists of 1 PDE. This function sets up the model parameters. 
% The function allows one input argument, which must be a struct. It is
% called inputStruct below, and can have the following fields (if not
% available then defaults are used):
%    inputStruct.params.selectICandBC  [1 or 2; default 1]
%                    1  ... initial condition: peak at top boundary,
%                           independent of x;  
%                           const Dirichlet BC at top, otherwise no flux
%                    2  ... initial condition: peak in middle,
%                           independent of x;  
%                           no flux BCs
%    inputStruct.params.yshift          [real >= 0; default 0]
%    inputStruct.params.rcells          [integer > 0; default 100]
%    inputStruct.params.selectReaction  ['none, 'logistic', 'genlogistic';
%                                        default 'none']
%    inputStruct.params.ReactionParams  [real vector with parameters;
%                                        default [] ]
%    inputStruct.grd.isAxiSymmetric     [true or false; default true]
%    inputStruct.tdr.tvec               [real vector with initial and 
%                                        output times; default [0:100]/100 ]
%
% Use ComparisonAxisymmetricAnd2D() to run the simulation.
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/Plasmin.m
%* Date created  : 2006, January 25
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     : 1.1 (AG, 2006 06 08)
%*                 - added inputStruct.params.rcells
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2006 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************


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
params.diffConst = [1];    % diffusion constants for the species

params.yshift = setValue('yshift', 0.0, inputStruct, 'params');
if (params.yshift < 0 )
  error('yshift must be non-negative!');
end

params.rcells =  setValue('rcells', 100, inputStruct, 'params');
if (params.rcells <= 0)
  error('rcells must be positive!');
end
if (params.rcells ~= round(params.rcells))
  error('rcells must be an integer!');
end

params.selectICandBC = setValue('selectICandBC', 1, inputStruct, 'params');

params.selectReaction = setValue('selectReaction', 'none', inputStruct, ...
                           'params'); 
params.ReactionParams = setValue('ReactionParams', [], inputStruct, 'params'); 

%
% params struct completed
%


%
% initialise struct tdr
%
tdr.size = 1; % number of species (PDEs)
% set output times
tdr.tvec = setValue('tvec', [0:100]/100, inputStruct, 'tdr'); 

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
%%% End of CODES %%%%


% describe reaction function ProbFReac.m 
tdr.FReac.depends = [ DependsU ];  % if only Zero entries then 
                                   % ProbFReac.m is never called
tdr.FReac.vec = true;  % true if ProbFReac.m can be called vectorized

% describe transport function ProbFTrans.m
tdr.FTrans.depends = [ Const ]; % if only Zero entries then 
                                % ProbFTrans.m is never called
%
% tdr struct completed
%


%
% define the spatial domain via its patches and store data in grd struct
%
grd.isAxiSymmetric = setValue('isAxiSymmetric', true, inputStruct, 'grd');

%%%ONE BIG PATCH%%%
%%% domain is y-shifted unit square
if (1)
  grd.nop = 1;             % Number Of Patches
  grd.ngb = [0 0 0 0];     % patchId of [left right bottom top] neighbour patch
                           % 0 means no neighbour in that direction
  patchId = 1;
  grd.nx(patchId) = 6;          % number of cells in x-direction 
  grd.ny(patchId) = params.rcells;     % number of cells in y-direction 
  grd.dx(patchId) = 1/grd.nx(patchId); % x-cell width of patch patchId
  grd.dy(patchId) = 1/grd.ny(patchId); % y-cell width of patch patchId
  grd.x0(patchId) = 0.0;           % x-coordinate of lower left corner of patch
  grd.y0(patchId) = params.yshift; % y-coordinate of lower left corner of patch
end
%%%End of ONE BIG PATCH%%%

           
%
% grd struct completed
%

return
% end of function


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
