function [a] = evalIntegral1D(g, mask, methodId)
% 
% This function evaluates the non-local term in 1D in all right
% cell end points for a grid function given in the column
% vector g. The result is returned as column vector a. If g is a
% matrix then the non-local term is evaluated for each column of g
% and the results are returned in a nmatrix a.
% For the evaluation two different algorithms can be
% employed and are selected by methodId 
%   methodId == 1 ... summation based
%   methodId == 2 ... FFT-based 
% Usually, methodId == 2 is the most efficient.
% The non-local term is encapsulated in the data structure mask,
% which is created by using functions setupIntegralRule1D_weights
% and setupIntegralRule1D_BCs.
%
% Alf Gerisch (gerisch@mathematik.tu-darmstadt.de)
% Version 1.1 of May 18, 2010
%  * addded support for nonperiodic bounday conditions.
% Version 1.0 of January 15, 2007
%


switch methodId
 case 1 
  % Non-FFT variant using weight vector w
  switch mask.BCs
   case 'pp'
    a=zeros(size(g));
    ge = [g(end-(mask.lm)+1:end,:); g; g(1:mask.lp,:)];
    w = mask.lm+mask.lp;
    % evaluate integral
    for i=1:size(g,1)
      a(i,:)=mask.weights*ge(i:(i+w),:);
    end
   otherwise
    error('evalIntegralRule1D:: unsupported value of BCs.');
  end
 case 2 
  % compute matrix-vector product via FFT
  switch mask.BCs
   case 'pp'
    a = ifft(repmat(mask.circFFT, 1, size(g,2)) .* fft(g));
   case 'zz'
    % extend g by zeros appropriately
    gext = [g; zeros(mask.N1ext-mask.N1,size(g,2))];
    aext = ifft(repmat(mask.circFFT, 1, size(g,2)) .* fft(gext));
    a    = aext(1:mask.N1+1,size(g,2));
   case 'vv'
    % here we assume that g is already extendend appropriately at
    % the left  (top   ; by mask.lm times size(g,2) entries) 
    % and right (bottom; by mask.lp times size(g,2) entries) end
    aext = ifft(repmat(mask.circFFT, 1, size(g,2)) .* fft(g));
    a    = aext(1:mask.N1,size(g,2));  
   case 'vz'
    % here we assume that g is already extendend appropriately at
    % the left  (top   ; by mask.lm+1 times size(g,2) entries) 
    % and extend further by zeros on the right (bottom)
    gext = [g; zeros(mask.N1ext-mask.N1-1-mask.lm,size(g,2))];
    aext = ifft(repmat(mask.circFFT, 1, size(g,2)) .* fft(gext));
    a    = aext(1:mask.N1+1,size(g,2));  
   case 'zv'
    % here we assume that g is already extendend appropriately at
    % the right (bottom; by mask.lp times size(g,2) entries) 
    % and extend further by zeros on the right (bottom)
    gext = [g; zeros(mask.lm,size(g,2))];
    aext = ifft(repmat(mask.circFFT, 1, size(g,2)) .* fft(gext));
    a    = aext(1:mask.N1,size(g,2));  
   otherwise
    error('evalIntegralRule1D:: unsupported value of BCs.');
  end
 otherwise
  error('evalIntegral1D::unknown methodId.');
end


return % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


