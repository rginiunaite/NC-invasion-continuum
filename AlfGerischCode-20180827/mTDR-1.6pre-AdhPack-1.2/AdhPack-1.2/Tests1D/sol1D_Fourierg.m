function a = sol1D_Fourierg(x, R, k, usesin, OmegaId)
% a = sol1D_Fourierg(x, R, k, usesin, OmegaId)
% 
% Function returns the exact value of the non-local term in 1D with
% sensing radius R and using Omega_i(r) with i=OmegaId evaluated at
% x for the function g(x) = sin(2*k*pi*x) if usesin==true and  g(x)
% = cos(2*k*pi*x) if usesin==false. 
%
% Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
% Version of 1.0 of January 15, 2007


switch OmegaId
 case 1
  if (usesin)
    a = (cos(2*k*pi*x) - cos(2*k*pi*(R + x)))/(4*k*pi*R) ...
	+ sin(k*pi*R)*sin(k*pi*(R - 2*x))/(2*k*pi*R);
    a = a/R;
  else
    a = -cos(k*pi*(R - 2*x))*sin(k*pi*R)/(2*k*pi*R) ...
	+ (-sin(2*k*pi*x) + sin(2*k*pi*(R + x)))/(4*k*pi*R);
    a = a/R;
  end
  
 case 2
  if (usesin)
    a = -(-4*k*pi*R*cos(2*k*pi*x) + sin(2*k*pi*(R - x)) + ...
	 sin(2*k*pi*(R + x)))/(4*(k*pi*R)^2);
    a = a/R;
  else
    a = -(-cos(2*k*pi*(R - x)) + cos(2*k*pi*(R + x)) ...
	  + 4*k*pi*R*sin(2*k*pi*x)) / (4*(k*pi*R)^2);
    a = a/R;
  end
 otherwise
  error('unknown OmegaId.')
end

return; % enf of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
