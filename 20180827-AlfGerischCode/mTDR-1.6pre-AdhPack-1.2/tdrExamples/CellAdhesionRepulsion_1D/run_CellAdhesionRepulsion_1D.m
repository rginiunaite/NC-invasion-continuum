function run_CellAdhesionRepulsion_1D(id)

if nargin==0
  id = 1;
end

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',20)

params.domainlength    = 1;
params.gridCells       = 400;

switch id
 case 1
  params.BCs = 'pp'; % periodic BCs
 otherwise
  error('unknown BC type');
end

inputStruct.params = params;
CellAdhesionRepulsion_1D(inputStruct);


