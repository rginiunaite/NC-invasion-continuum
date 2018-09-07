function [t, y, stats, hpred] = rowmapStepOn(RHSFun, tspan, y0, options)
% [t, y, stats, hpred] = rowmapStepOn(RHSFun, tspan, y0, options)
% 
% rowmapStepOn calls ROWMAP repeatedly to step onto the values
% provided in tspan exactly (avoiding interpolation). This is
% useful for problems with discontinuous right-hand side function
% where the discontinuity is at predefined time points.
%
% RHSFun, tspan, y0 are mandatory input arguments, options is
% optional. 
%
% The option ReturnMode is ignored with this call. If t and y are
% requested return values then t=tspan and the corresponding
% solution approximations are in the rows of y (this corresponds to
% ReturnMode==2); if t and y are not requested then the function
% behaves like for ReturnMode==0. The user can use a
% output and/or post-step function to store other solution
% approximations. Note that these functions, if provided,
% are called with flag=='init' or flag=='done' repeatedly (once for
% each tspan(1:end-1)). The return arguments stats and hpred are
% also optional.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File    : rowmapStepOn.m
% Version : 25 January 2009 (Alf Gerisch, University of Halle)
%
% Documentation of the ROWMAP MEX Interface is maintained at
% http://sim.mathematik.uni-halle.de/~gerisch/r/rowmapmex.html
%
% Please send comments and bug reports to 
%     alf.gerisch@mathematik.uni-halle.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin<4) % we have no options structure
  options.ReturnMode = 1;
else % we have an options structure
  if isfield(options, 'ReturnMode') 
    if ((options.ReturnMode ~= 2) && (nargout >= 2))
      disp(['rowmapStepOn:: Ignoring provided option value ReturnMode = '...
	    num2str(options.ReturnMode) '.'])
    end
  end
  options.ReturnMode = 1;
end

if (nargout >= 2)
  t=tspan;
  y=zeros(length(t),length(y0));
  y(1,:) = y0';
end
y0_i = y0; % initial value for next call to rowmap

% initial time step size 
if isfield(options, 'InitialStep')
    hpred = options.InitialStep;
end

% initialise stats vector
stats = zeros(1,5);

for ii=2:length(tspan)
  tspan_i = tspan([ii-1 ii]);
  options.InitialStep = hpred; % set initial time step size
  [ti, yi, statsi, hpred] = rowmap(RHSFun, tspan_i, y0_i, options);
  if (nargout >= 2) % store computed solution
    t(ii) = ti;
    y(ii,:) = yi';
  end
  y0_i = yi'; % initial value for next call to rowmap
  % update stats vector
  stats(1) = statsi(1);
  stats(2) = stats(2) + statsi(2);
  stats(3) = stats(3) + statsi(3);
  stats(4) = stats(4) + statsi(4);
  stats(5) = stats(5) + statsi(5);
end