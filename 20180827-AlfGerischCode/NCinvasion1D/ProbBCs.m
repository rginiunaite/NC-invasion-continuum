%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [BCtype,BCval] = ProbBCs(patchId, boundary, t, xvals, yvals, params)
%
% Input arguments
%   patchId       ... [integer, >=1] characterising the patch for which 
%                     boundary data is requested.
%   boundary      ... [one of 1,2,3,4] characterising the boundary (left,
%                     right, bottom, top) of patch patchId for which
%                     boundary data is requested.
%   t             ... current time
%   xvals,yvals   ... [real vectors of same length] (xvals(i), yvals(i)) are 
%                     the coordinates for which boundary data is requested.
%   params        ... [struct] as provided by ProbGetParams().
% Output arguments
%   BCtype, BCval ... [integer and real matrix of size length(xvals) times 
%                     no of PDEs in system] The (i,j)th entry of BCtype 
%                     encodes the type of boundary condition in (xvals(i), 
%                     yvals(i)) of equation j. Boundary condition types are
%                     None      = 0;
%                     ZeroFlux  = 1;
%                     Dirichlet = 2;
%                     If BCtype(i,j) == Dirichlet then BCval(i,j) contains the
%                     corresponding Dirichlet value, otherwise BCval(i,j)
%                     is unused.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BCtype,BCval] = ProbBCs(patchId, boundary, t, yvals, params)


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


if boundary == left
    % left BCs
    BCtype = [ZeroFlux Dirichlet ZeroFlux];
    BCval = nan(size(BCtype));
    BCval(1,2) = 1;
elseif boundary == right
    % right BCs
    BCtype = [ZeroFlux ZeroFlux ZeroFlux];
    BCval = nan(size(BCtype));
else
    error('what???')
end

return;
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

