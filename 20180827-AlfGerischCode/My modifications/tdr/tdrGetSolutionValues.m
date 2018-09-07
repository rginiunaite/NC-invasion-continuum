function sol = tdrGetSolutionValues(y, patchId, rows, cols, spec)

global TDRP;

if (TDRP.grd.is1D)
  error('tdrGetSolutionValues::not implemented for 1D space.');
end

% Determine ystart such that y(ystart) corresponds to the first
% solution value of patch patchId in y
ystart = TDRP.grd.ps(patchId);

% compute indices into y for rows and cols for species 1
colstarts = [ystart + (cols(:)-1)*TDRP.grd.ny(patchId)*TDRP.tdr.size]';
rowshift = (rows(:)-1)*TDRP.tdr.size;
ind = ones(length(rows),1)*colstarts+rowshift*ones(1,length(cols));
ind = ind(:);

% fill return value sol
if (size(spec) == [1,1])
  sol = reshape(y(ind+spec-1), length(rows), length(cols));
else
  sol = NaN(length(rows), length(cols), length(spec));
  for ss=1:length(spec)
    sol(:,:,ss) = reshape(y(ind+spec(ss)-1), length(rows), length(cols));
  end
end  

return;
% end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

