function run_CellAdhesion_2D_1species(id)

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)

params.domainlength    = 1;
params.gridCells       = 150;

switch id
 case 1
  params.BCs = 'pppp';  % periodic in both directions
 case 2
  params.BCs = 'zzzz';  % zero-flux on all of boundary
 otherwise
  error('unknown BC type');
end

inputStruct.params = params;
CellAdhesion_2D_1species(inputStruct);


