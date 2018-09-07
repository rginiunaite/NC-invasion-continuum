function nonlocalBCs = setupIntegralRule1D_BCs(nonlocal, N1, BCs)
% function nonlocalBCs = setupIntegralRule1D_BCs(nonlocal, BCs)
%
% Add boundary condition information to nonlocal structure
% (computed from setupIntegralRule1D_weights.m). 
% N1          ... number of grid points (cells) in x1 direction
% 
% The following fields are added to nonlocal
%   BCs         ... the boundary condition to use
%   N1          ... copied from input
%   N1ext       ... dimension of the circulant matrix circ
%   circ        ... first column of the circulant matrix circ
%   circFFT     ... FFT of circ
%
% Valid values for BCs are 'pp' ... periodic
%                          'zz' ... zero left and right
%                          'zv' ... zero left and values right
%                          'vz' ... values left and zero right
%                          'vv' ... values left and right
%
% Alf Gerisch (gerisch@mathematik.tu-darmstadt.de)
% Version 1.0 of May 18, 2010
%   * initial version
%

% grab some values out of nonlocal
lm = nonlocal.lm;
lp = nonlocal.lp;
weights = nonlocal.weights;

switch BCs
 case 'pp'
  % the matrix is a banded circulant matrix
  circ = [weights(lm+1:-1:1) ...
          zeros(1, N1-length(weights)) ...
          weights(end:-1:lm+2)]';
  N1ext = N1;
 case 'zz'
  % the matrix is a banded Toeplitz matrix
  Trows = N1+1;
  Tcols = N1;
  Tlp   = lp;
  Tlm   = lm;
  N1ext = max(Trows+Tlm, Tcols+Tlp);
  circ = [weights(lm+2:-1:1) ... 
          zeros(1, N1ext-lm-lp-1) ...
          weights(end:-1:lm+3)]';
 case 'vv'
  % the matrix is an upper triangular, Toeplitz matrix
  N1ext = N1+lp+lm;
  circ = [weights(1) ...
	  zeros(1, N1ext-lm-lp-1) ...
	  weights(end:-1:2)]';
 case 'vz'
  % the matrix is an upper triangular, Toeplitz matrix
  Trows = N1+1;
  Tcols = N1+1+lm;
  Tlp   = lp;
  Tlm   = lm;
  N1ext = max(Trows+Tlm, Tcols+Tlp);
  circ = [weights(1) ...
	  zeros(1, N1ext-lm-lp-1) ...
	  weights(end:-1:2)]';
 case 'zv'
  % the matrix is a Toeplitz matrix
  N1ext = N1+lm+lp;
  circ = [weights(lm+1:-1:1) ...
	  zeros(1, N1ext-length(weights)) ...
	  weights(end:-1:lm+2)]';
  length(circ)
 otherwise
  error('setupIntegralRule1D_BCs:: unsupported value of BCs provided.')
end


%{
aa=circulant(circ)*[ones(N1,1);zeros(N1ext-N1,1)];
if strcmp(BCs, 'pp')
  aa=aa(1:N1);
else
  aa=aa(1:N1+1);
end
aa([1:20:200])
aa([400:200:end-200])
aa([end-200:20:end])

error('brrr');
%}

nonlocalBCs         = nonlocal;
nonlocalBCs.N1      = N1;
nonlocalBCs.BCs     = BCs;
nonlocalBCs.N1ext   = N1ext;
nonlocalBCs.circ    = circ;
nonlocalBCs.circFFT = fft(nonlocalBCs.circ);

return % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

