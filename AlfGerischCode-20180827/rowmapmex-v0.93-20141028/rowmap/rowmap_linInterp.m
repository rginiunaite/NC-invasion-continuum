function yi = rowmap_linInterp(told, yold, tnew, ynew, t)
% function yi = rowmap_linInterp(told, yold, tnew, ynew, t)
%
% This function provides a Matlab implementation for interpolating
% ROWMAP solution approximations. It can be used in a post-step or
% output function called from rowmap. 
%
% This function implements linear interpolation.
%
% It should be called only for t values satisfying told <= t <=
% tnew. yi is the interpolated value at t using the data yold at
% told and ynew at tnew.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : rowmap_linInterp.m
% Version : 25 January 2009 (Alf Gerisch, University of Halle)
%
% Documentation of the ROWMAP MEX Interface is maintained at
% http://sim.mathematik.uni-halle.de/~gerisch/r/rowmapmex.html
%
% Please send comments and bug reports to 
%     alf.gerisch@mathematik.uni-halle.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fac = (t - told)/(tnew-told);
yi = (1-fac)*yold + fac*ynew;
return;

