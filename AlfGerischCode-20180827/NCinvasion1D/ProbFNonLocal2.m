function [P] = ProbFNonLocal2(i,params,varargin)


switch i
 case params.eq.n1
  n1 = varargin{1}(:,params.eq.n1);
  n2 = varargin{1}(:,params.eq.n2);
  P = ones(size(n1));
 case params.eq.n2
  n1 = varargin{1}(:,params.eq.n1);
  n2 = varargin{1}(:,params.eq.n2);
  P = ones(size(n1));
 otherwise
  error('Subroutine should not have been called.');
end
return


% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
