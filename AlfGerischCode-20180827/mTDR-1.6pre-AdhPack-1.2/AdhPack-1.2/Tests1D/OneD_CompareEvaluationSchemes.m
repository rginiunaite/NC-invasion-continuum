function OneD_CompareEvaluationSchemes(runOrPlot)
%{
 diary('OneD_CompareEvaluationSchemes.diary');
 OneD_CompareEvaluationSchemes('run');
 OneD_CompareEvaluationSchemes('plot')
 printeps('OneD_CompareEvaluationSchemes.eps');
 diary off
%} 






oldpath = path;
addpath('..');

format compact
format short e

d=[];
if strcmp(runOrPlot, 'run')
  Rvals  = [0.05 0.1 0.2 0.4 0.8]
  Nvals = round(1./logspace(log10(1/5),log10(1/10005),30))
  OmegaId = 2;
  RuleId  = 2;
  NR = 3000;
  
  d.timings = [];
  d.Nvals   = Nvals;
  d.Rvals   = Rvals;
  d.OmegaId = OmegaId;
  d.RuleId  = RuleId;
  d.NR      = NR;
  d.c1=2;

  for iN = 1:length(Nvals)
    N = Nvals(iN);
    h = 1/N;
    N1=d.c1*N;
    x1=([0:(N1-1)]'+0.5)*h;
    for iR = 1:length(Rvals)
      R = Rvals(iR);
      disp(['(N,R) = ' num2str([N,R])]);
      % setup integration rule
      mask = setupIntegralRule1D(R, h, N1, OmegaId, RuleId, NR);
      % create 400 random vectors
      g = rand(N1,400);
      % evaluate non-local term and time schemes
      for evalId = [1 2] % 1=summation, 2=FFT
	tstart=cputime; 
	res{evalId} = evalIntegral1D(g, mask, evalId);
	elapsed = cputime-tstart;
	d.timings{iN}{iR}{evalId} = elapsed;
      end % for evalId
      % compute the difference between the results of the two
      % evaluation schemes
      err = norm(res{1}-res{2},inf);
      disp(['Difference between results of evaluation schemes = ' num2str(err)]);
      if (err>1e-12)
        warning('difference in evaluation schemes larger than 1e-12.')
      end
    end % for iR = 1:length(Rvals)
  end %for iN = 1:length(Nvals)

  
  save('OneD_CompareEvaluationSchemes.mat','d');
  rmpath('..')
  return
end


if strcmp(runOrPlot, 'plot')

  load('OneD_CompareEvaluationSchemes.mat')

  %set(0,'DefaultAxesFontSize',20)
  %set(0,'DefaultTextFontSize',20)

  set(0,'DefaultAxesFontSize',15)
  set(0,'DefaultTextFontSize',15)
  set(0,'DefaultLineLineWidth',1)
  set(0,'DefaultAxesLineWidth',1)
  set(0,'DefaultTextInterpreter','latex');
  figure(1)
  
  timesum = zeros(length(d.timings),1);
  timefft = timesum;
  for iR = 1:length(d.Rvals);
    for iN = 1:length(d.timings)
      timesum(iN) = d.timings{iN}{iR}{1};
      timefft(iN) = d.timings{iN}{iR}{2};
    end
    %[timesum timefft]
    lh = loglog(d.Nvals,timesum, 'g-s', ...
		d.Nvals, timefft, 'r-x');
    hold on
   % pause
  end  
    
  c0 = 2*d.Rvals(3);
  loglog(d.Nvals, 1.5e-5*d.c1*c0*(d.Nvals.^2),'--', ...
	 'color', 0.7*[0.5 1 0.5],'linewidth',2)
  loglog(d.Nvals, 3e-6*2*d.c1*d.Nvals.*log2(d.Nvals),'--', ...
	 'color', 0.7*[1 0.5 0.5],'linewidth',2);
  axis([10, 11000, 5e-3 800])
  xlabel('$N$'), 
  ylabel('CPU time (seconds)')

  [LEGH,OBJH,OUTH,OUTM] = ...
      legend(lh, {'summation', 'FFT'}, 'location', 'northwest');
  for hhi = 1:length(OBJH)
    hh = OBJH(hhi);
    if strcmp(get(hh,'type'), 'text')
      set(hh,'interpreter','latex');
    end
  end

  rmpath('..')
  return
end





