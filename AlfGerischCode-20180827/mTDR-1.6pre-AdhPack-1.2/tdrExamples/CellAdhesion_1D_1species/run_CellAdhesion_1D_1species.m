function run_CellAdhesion_1D_1species(id)

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)

params.domainlength    = 1;
params.gridCells       = 400;

switch id
 case 1
  params.BCs = 'pp'; % periodic BCs
 case 2
  params.BCs = 'vz'; % symmetry at left and zero-flux at right end
 case 3
  params.BCs = 'zz'; % zero flux BCs on both ends
 otherwise
  error('unknown BC type');
end

inputStruct.params = params;
CellAdhesion_1D_1species(inputStruct);


