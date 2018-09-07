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
%   varargin      ... 
%                     {1} contains time
%                     {2} contains the results for the variables
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

t = varargin{1};
Lt2 = domaingrowth(t,params)^2

if ((i == params.eq.u) & (j == params.eq.u))
  pij = params.u_D/Lt2;
elseif ((i == params.eq.v) & (j == params.eq.v))
  pij = params.v_D/Lt2;
elseif ((i == params.eq.c) & (j == params.eq.c))
  pij = params.c_D/Lt2;
elseif ((i == params.eq.u) & (j == params.eq.c))
  u=varargin{2}(:,params.eq.u);
  v=varargin{2}(:,params.eq.v);
  c=varargin{2}(:,params.eq.c);
  switch params.selectChemotaxisLeader
      case 'receptor'
          pij = params.chemotaxisLeader_k*params.chemotaxisLeader_beta ...
              ./ ((c+params.chemotaxisLeader_beta).^2 * Lt2) ...
             .* max(0,1-u-v);
           
          % add volume filling here: .* (1-u-v)
      case 'simple'
          pij = params.chemotaxisLeader_k ...
              .* max(0,1-u-v);
      otherwise
          error('undefined params.selectChemotaxisLeader');
  end
elseif ((i == params.eq.v) & (j == params.eq.u))
  u=varargin{2}(:,params.eq.u);
  v=varargin{2}(:,params.eq.v);
  switch params.selectChemotaxisFollower
      case 'receptor'
          pij = params.chemotaxisFollower_k*params.chemotaxisFollower_beta ...
              ./ ((u+params.chemotaxisFollower_beta).^2 * Lt2) ...
              .* max(0,1-u-v);
          % add volume filling here: .* (1-u-v)
      case 'simple'
          pij = params.chemotaxisFollower_k ...
              .* max(0,1-u-v);      
      otherwise
          error('undefined params.selectChemotaxisFollower');
  end
else
  error(['ProbFTrans(' num2str(i) ',' num2str(j) ...
	 ') should not have been called.']);
end

return
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
