function C = circulant(c)
% function C = circulant(c)
%
% Generate the full circulant matrix C defined by its first column c.
%
% A. Gerisch (gerisch@mathematik.tu-darmstadt.de, 2010-05-13)
%

% make sure c is a column
c=c(:);
% get dimension
n = length(c);
% preallocate C
C = NaN(n,n);
% fill C
for ii=1:n
  C(:,ii)=circshift(c,ii-1);
end

return
%end of function
