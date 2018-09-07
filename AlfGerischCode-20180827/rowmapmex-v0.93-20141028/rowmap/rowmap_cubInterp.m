function yi = rowmap_cubInterp(told, yold, fold, tnew, ynew, fnew, t)
% function yi = rowmap_cubInterp(told, yold, fold, tnew, ynew, fnew, t)
%
% This function provides a Matlab implementation for interpolating
% ROWMAP solution approximations. It can be used in a post-step or
% output function called from ROWMAP. 
%
% This function implements cubic (Hermite) interpolation.
%
% It should be called only for t values satisfying told <= t <=
% tnew. yi is the interpolated value at t using the data yold, fold at
% told and ynew, fnew at tnew.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : rowmap_cubInterp.m
% Version : 25 January 2009 (Alf Gerisch, University of Halle)
%
% Documentation of the ROWMAP MEX Interface is maintained at
% http://sim.mathematik.uni-halle.de/~gerisch/r/rowmapmex.html
%
% Please send comments and bug reports to 
%     alf.gerisch@mathematik.uni-halle.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hs      = tnew-told;
theta   = (t-told)/hs;

g1 = theta*theta*(3-2*theta);
g2 = hs*theta*(theta-1)^2;
g3 = hs*theta*theta*(theta-1);
yi = (1-g1)*yold+g1*ynew+g2*fold+g3*fnew;

return;

