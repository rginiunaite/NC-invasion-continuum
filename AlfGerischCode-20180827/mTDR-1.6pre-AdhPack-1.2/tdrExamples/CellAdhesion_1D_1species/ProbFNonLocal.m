function [G] = ProbFNonLocal(i,params,varargin)


switch i
 case params.eq.c1
  % varargin{1}(:,params.eq.c1) is the c1 density in all grid cells
  G = (params.c_Sc1c1*varargin{1}(:,params.eq.c1)) ...
      .* max(0,1-varargin{1}(:,params.eq.c1));
 otherwise
  error('Subroutine should not have been called.');
end
return


% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
