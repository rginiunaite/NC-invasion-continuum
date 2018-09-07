function [oldpath] = tdrSetup()
% tdrSetup()
%
% This function sets the paths for the mTDR system.
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrSetup.m
%* Date created  : 2008, September 8
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0 (new in Revison 1.3 of mTDR)
%* Revisions     :
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

oldpath = path;

% the primary mTDR subdirectory is in the path already when this
% function is called, so we try to get it here
TDRbasedir = which('tdrSetup');
TDRbasedir = TDRbasedir(1:(end-length('tdrSetup.m')));
disp(['tdrSetup:: using mTDR base directory ' TDRbasedir]);

% add absolute base directory and subdirectories to path
path(TDRbasedir, path)
path([TDRbasedir '/tdr'], path);
path([TDRbasedir '/tdrUtil'], path);
path([TDRbasedir '/AdhPack-1.2'], path);

return
%end of function
