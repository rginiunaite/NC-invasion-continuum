%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function pij = ProbFTrans(i, j, params, varargin)
%
% Input arguments
%   i             ... variable who's transport phenomena are being
%                     described
%   j             ... variable who's concentration/density is determining
%                     the transport of variable i
%   params        ... [struct] as provided by ProbGetParams().
%   varargin      ... {1} contains the results for the variables
%                          -> varargin{1}(:,:,c_m) is a matrix containing
%                          the concentration of MSC in all x and y 
%                     {2} is the patchId 
%                     {3} is either x or y
% Output arguments: None
% 
% Return function values of parameter function p_{i,j}(...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pij, terms] = ProbFTrans(i, j, params, varargin)


if ((i == params.eq.c1) & (j == params.eq.c1))
  pij = params.c1_D;
else
  error(['ProbFTrans(' num2str(i) ',' num2str(j) ...
	 ') should not have been called.']);
end

return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
