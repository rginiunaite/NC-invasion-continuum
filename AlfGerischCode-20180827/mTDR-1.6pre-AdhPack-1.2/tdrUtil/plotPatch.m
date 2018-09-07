function [patchHdl] = plotPatch(x0,dx,nx,y0,dy,ny,zvals,params_)
%
% function [patchHdl] = plotPatch(x0,dx,nx,y0,dy,ny,zvals,params)
%
% This function plots a 3D patch of data given as cell averages on a
% uniform grid of cells in a rectangular area.
%
% Input:
%    x0, y0   ... lower left corner of patch
%    dx, dy   ... x- and y-grid width (dimension of the cells)
%    nx, ny   ... number of grid cells in x- and y-direction
%    zvals    ... (ny,nx) matrix of average z-values in all grid cells
%    params   ... struct of additional parameters (optional, defaults
%                 below)
%
% Output:
%    patchHdl ... Patch handle or matrix of patch handles created
%                 (depending on selected method).
%
% Additional Parameters:
%    params.methodId  ... selects the plot method (default is 1)  
%       methodId=1: Fast Method: compute average of zval in the corners of
%                   the grid cells by using the 4 adjacent grid cells of
%                   such a corner (accordingly for the boundary of the
%                   patch) and plot one big patch using this data (hence
%                   only one patch handle is returned). This method might
%                   result in some parts of the patch not displaying
%                   properly in view(2) mode. 
%       methodId=2: Slower method: generate for each grid cell a patch
%                   according to the corresponding height in zval. This
%                   results in a mtarix of patch handles, which are
%                   returned. This method is much slower but the result in
%                   view(2) mode is better (no parts of the patch are
%                   missing). 
%
%-----!-----------------------------------------------------------------
%     AUTHOR:   ALF GERISCH (GERISCH@MATHEMATIK.UNI-HALLE.DE)
%               FB MATHEMATIK / INFORMATIK
%               MARTIN-LUTHER-UNIVERSIT"AT HALLE-WITTENBERG
%               POSTFACH, 06099 HALLE (SAALE), GERMANY
%-----!-----------------------------------------------------------------

% check argument list
if (nargin<7) 
  error('Need at least 7 input arguments'); 
end

% copy additional parameter if supplied) to struct params
if (nargin>7) 
  params = params_; 
else 
  params.empty=true; % dummy entry
end

% check additional parameter

% params.methodId
fieldName='methodId';
allowedValues=[1 2];
default = 1;
if (isfield(params, fieldName)) 
  if (min(abs(params.(fieldName)-allowedValues))>0)
    error(['Only values "' num2str(allowedValues) '" are allowed for ' ...
       fieldName '. You supplied ' num2str(params.(fieldName)) '.']);
  end
else
  params.(fieldName) = default;  
end


if (params.methodId==1) 

  %%%%%%%%%%%
  % Method 1%
  %%%%%%%%%%%
  x=[x0 x0+dx/2:dx:(x0+(nx-1)*dx) x0+nx*dx];
  y=[y0 y0+dy/2:dy:(y0+(ny-1)*dy) y0+ny*dy];
  [XX,YY]=meshgrid(x,y);
  
  ZZ=zeros(size(XX));
  ZZ(2:(end-1),2:(end-1))=(zvals(1:end-1, 1:end-1) + ...
               zvals(2:end,   1:end-1) + ...
               zvals(1:end-1, 2:end) + ...
               zvals(2:end,   2:end))/4;
  ZZ(1,1) = zvals(1,1);
  ZZ(1,2:end-1) = (zvals(1,1:end-1)+zvals(1,2:end))/2;
  ZZ(1,end) = zvals(1,end);
  ZZ(2:end-1,end) = (zvals(1:end-1,end)+zvals(2:end,end))/2;
  ZZ(end,end) = zvals(end,end);
  ZZ(end,2:end-1) = (zvals(end,1:end-1)+zvals(end,2:end))/2;
  ZZ(end,1) = zvals(end,1);
  ZZ(2:end-1,1) =(zvals(1:end-1,1)+zvals(2:end,1))/2;
  
  patchHdl=surf(XX,YY,ZZ);
  shading interp
else

 %%%%%%%%%%%
 % Method 2%
 %%%%%%%%%%%
 x=[x0:dx:(x0+nx*dx)];
 y=[y0:dy:(y0+ny*dy)];
 
 patchHdl=zeros(size(zvals));
 for i=1:nx
   for j=1:ny
     patchHdl(j,i)=patch([x(i) x(i+1) x(i+1) x(i)], ...
            [y(j) y(j) y(j+1) y(j+1)], ...
            zvals(j,i) * [ 1 1 1 1],...
            zvals(j,i), 'LineStyle', 'none');
   end
 end
 
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $Header: /home/gerisch/CVS/Bone/BONE-EXAMPLE/MATLAB/plotPatch.m,v 1.2 2004/05/25 11:53:58 gerisch Exp $
%
% $Log: plotPatch.m,v $
% Revision 1.2  2004/05/25 11:53:58  gerisch
%     * minor changes, syntax corrections etc.
%
% Revision 1.1  2004/05/25 11:08:52  gerisch
%     * function [patchHdl] = plotPatch(x0,dx,nx,y0,dy,ny,zvals,params)
%     This function plots a 3D patch of data giving as cell averages
%     on a uniform grid of cells in an rectangular area.
%
%
