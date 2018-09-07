function [A11, A22] = evalIntegral2D_old(G, mask, methodId)
% [A11, A22] = evalIntegral2D(G, mask, methodId)
%
% This function evaluates the first component of the non-local term
% in 2D in the centre of all right cell faces (returned as matrix
% A11) and its second component in the centre of all upper cell
% faces  (returned as matrix A22) for a periodic (in both
% directions) grid function G of dimension mask.x2Dir.N2 times
% mask.x2Dir.N1 (alternatively, G can be a vector of length
% mask.x2Dir.N2*mask.x2Dir.N1, which is formed columenwise from
% that matrix). For the evaluation three different algorithms can
% be employed and are selected by methodId 
%   methodId == 1 ... summation based
%   methodId == 2 ... FFT-based I
%   methodId == 3 ... FFT-based II
% Usually, methodId == 3 is the most efficient.
% The non-local term is encapsulated in the data structure mask,
% which must contain two subfields x1Dir and x2Dir corresponding to
% the right and upper cell faces, respectively. Using function
% setupIntegralRule2D, this structure is set up as
%   mask.x1Dir = setupIntegralRule2D('x1',...);
%   mask.x2Dir = setupIntegralRule2D('x2',...);
% These calls also prepare data within the returned structures
% which is required for the FFT-based evaluation schemes.
%
% Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
% Version 1.0 of January 15, 2007
% Version 1.1 of Sep 808, 2008
%  * allowed for G to be a matrix and not as before just a vector.
%

if (size(G) ~= [mask.x2Dir.N2, mask.x2Dir.N1])
  % G must be a vector, so reshape to matrix
  G = reshape(G, mask.x2Dir.N2, mask.x2Dir.N1);
end

switch methodId
 case 1  
  % Non-FFT Variant
  % extend for easy computation of integral (periodic BCs)
  km = max(mask.x1Dir.km,mask.x2Dir.km);
  kp = max(mask.x1Dir.kp,mask.x2Dir.kp);
  lm = max(mask.x1Dir.lm,mask.x2Dir.lm);
  lp = max(mask.x1Dir.lp,mask.x2Dir.lp);
  Ge = [
      G(end-km+1:end,end-lm+1:end) G(end-km+1:end,:) G(end-km+1:end,1:lp)
      G(        :   ,end-lm+1:end) G                 G(        :   ,1:lp)
      G(       1:kp ,end-lm+1:end) G(        1:kp,:) G(       1:kp ,1:lp)
       ];
  % evaluate integral
  A11=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
  A22=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
  yindices1 = km +[-mask.x1Dir.km:mask.x1Dir.kp];
  yindices2 = km +[-mask.x2Dir.km:mask.x2Dir.kp];
  xindices1 = lm +[-mask.x1Dir.lm:mask.x1Dir.lp];
  xindices2 = lm +[-mask.x2Dir.lm:mask.x2Dir.lp];
  for i=1:mask.x2Dir.N1
    for j=1:mask.x2Dir.N2
      tmp = mask.x1Dir.Wx1.*Ge(j+yindices1,i+xindices1);
      A11(j,i)= sum(tmp(:));
      tmp = mask.x2Dir.Wx2.*Ge(j+yindices2,i+xindices2);
      A22(j,i)= sum(tmp(:));
    end
  end

 case 2  
  % FFT Variant I:  
  % compute fft of G
  gFFT = fft(G);
  % construct A11 and A22 matrices
  A11=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
  A22=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
  for i=1:mask.x2Dir.N1
    A11(:,i) = sum(gFFT.*mask.x1Dir.Vx1FFT.viFFT(:,[i:-1:1 end:-1:i+1]), 2);
    A22(:,i) = sum(gFFT.*mask.x2Dir.Vx2FFT.viFFT(:,[i:-1:1 end:-1:i+1]), 2);
  end
  A11 = ifft(A11);
  A22 = ifft(A22);
  
 case 3  
  % FFT Variant II:
  % compute fft of G and a second fft to get the u_j vectors as columns in u
  gFFT = fft(G);
  u = gFFT(:);
  u = fft(reshape(u(mask.x1Dir.Vx1FFT.permv), mask.x2Dir.N1, mask.x2Dir.N2));  
  
  % multiply component wise with appropriate deltajFFT and perform ifft
  s1 = ifft(u .* mask.x1Dir.Vx1FFT.deltajFFT);
  s2 = ifft(u .* mask.x2Dir.Vx2FFT.deltajFFT);
  
  % make columns, permute, make matrix and perform ifft 
  s1 = s1(:);
  s2 = s2(:);
  A11 = ifft(reshape(s1(mask.x1Dir.Vx1FFT.permv), mask.x2Dir.N2, mask.x2Dir.N1));
  A22 = ifft(reshape(s2(mask.x1Dir.Vx1FFT.permv), mask.x2Dir.N2, mask.x2Dir.N1));
 
 otherwise
  error('evalIntegral2D::Matrix-Vector product methodId non-existent.');
end

return; % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  

