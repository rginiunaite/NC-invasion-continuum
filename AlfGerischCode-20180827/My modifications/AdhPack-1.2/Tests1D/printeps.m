function printeps(epsname)
% function printeps(epsname)
%
% print current figure to eps-file epsname (the name must be with
% extension) and add bounding box to the beginning. The printing is done in
% painters-mode.  
%
% This is a work-around to a current bug of Matlab (Version 7.3.0.298
% (R2006b), linux version) when ps-fonst are included in the eps-file. The
% bug prevents a successful execution of epstopdf to the generated
% eps-file. Works only on a Unix system!
%
% Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
% Version of 1.0 of June 15, 2007
%

printmode = 'painters';
disp(['printeps:: printing in ' printmode '-mode to ' epsname]);
print('-depsc', '-painters', epsname);
disp(['printeps:: writing bounding box info in first lines of file.']);
! echo '%!PS-Adobe-2.0 EPSF-1.2' > xx.eps
eval(['! echo $(grep ''%%BoundingBox'' ' epsname ' ) >> xx.eps'])
eval(['! cat ' epsname ' >> xx.eps']) 
eval(['!mv xx.eps ' epsname])
disp(['printeps:: done.']);

return;
