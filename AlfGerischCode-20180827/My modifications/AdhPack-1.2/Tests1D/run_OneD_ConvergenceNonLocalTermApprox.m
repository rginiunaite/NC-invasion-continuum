diary('run_OneD_ConvergenceNonLocalTermApprox.diary')
  
R       = 0.1;
  
for NR = [250 500 2000]
for OmegaId=[1,2]
  
  sink    = [2];
  cosk    = [];
  fname = OneD_ConvergenceNonLocalTermApprox(R, OmegaId, sink, cosk, NR);
  figure(1);printeps([fname '_RuleId1.eps'])
  figure(10);printeps([fname '_RuleIdAll.eps'])
 
  close all
  
  sink    = [];
  cosk    = [2];
  fname = OneD_ConvergenceNonLocalTermApprox(R, OmegaId, sink, cosk, NR);
  %figure(1);printeps([fname '_RuleId1.eps'])
  figure(10);printeps([fname '_RuleIdAll.eps'])
 
  close all
  
  sink    = [8];
  cosk    = [];
  fname = OneD_ConvergenceNonLocalTermApprox(R, OmegaId, sink, cosk, NR);
  %figure(1);printeps([fname '_RuleId1.eps'])
  figure(10);printeps([fname '_RuleIdAll.eps'])
 
  close all
    
  sink    = [2];
  cosk    = [2];
  fname = OneD_ConvergenceNonLocalTermApprox(R, OmegaId, sink, cosk, NR);
  %figure(1);printeps([fname '_RuleId1.eps'])
  figure(10);printeps([fname '_RuleIdAll.eps'])
 
  close all

end %for OmegaId=[...]
end %for NR = [...]

diary off
return
  
 