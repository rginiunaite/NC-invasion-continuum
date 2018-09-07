function TwoD_CompareEvaluationSchemes(runOrPlot)
%{
 diary('TwoD_CompareEvaluationSchemes.diary');
 TwoD_CompareEvaluationSchemes('run');
 TwoD_CompareEvaluationSchemes('plot')
 printeps('TwoD_CompareEvaluationSchemes.eps');
 diary off
%} 

format compact
format short e



oldpath = path;
addpath('..');

if strcmp(runOrPlot, 'run')

  Rvals  = [0.05 0.1 0.2 0.4 0.8]
  Nvals = round(1./logspace(log10(1/5),log10(1/800),20))
  %Nvals = round(1./logspace(log10(1/5),log10(1/400),10))
  OmegaId = 2  % 1 or 2
  RuleId  = 1  % 1 or 2
  
  d.timings = [];
  d.Nvals   = Nvals;
  d.Rvals   = Rvals;
  d.OmegaId = OmegaId;
  d.RuleId  = RuleId;
  d.c1=2;
  d.c2=2;

  for iN = 1:length(Nvals)
    N = Nvals(iN);
    h = 1/N;
    N1=d.c1*N;
    N2=d.c2*N;
    x1=([0:(N1-1)]'+0.5)*h;
    x2=([0:(N2-1)]'+0.5)*h;
    for iR = 1:length(Rvals)
      R = Rvals(iR);
      disp(['(N,R) = ' num2str([N,R])]);
      % setup integration rule
      switch OmegaId
       case 1
	nl.OmegaStr = '@(r)(1/(pi*R^2))';
       case 2
	nl.OmegaStr = '@(r)(3/(pi*R^2)*(1-r/R))'; 
       otherwise
	error('unknown case.')
      end
      nl.x1Dir = setupIntegralRule2D('x1', R, h, N1, N2, eval(nl.OmegaStr), ...
					   RuleId, 100, +inf);
      nl.x2Dir = setupIntegralRule2D('x2', R, h, N1, N2, eval(nl.OmegaStr), ...
					   RuleId, 100, +inf);
      disp('... rule setup complete')
      
      % create 10 random vectors
      g = rand(N1*N2,10);
      
      % select evaluation schemes depending on problem size
      evalIds = [1 2 3];
      if (N>200)
	evalIds = [2 3];
      end
      if (N>500)
	evalIds = [3];
      end
      
      % evaluate non-local term and time schemes
      for evalId = evalIds % 1=summation, 2=FFTi 3=FFTii
	tstart=cputime; 
	for j=1:size(g,2)
	  [res1, res2] = evalIntegral2D(g(:,j), nl, evalId);
	end
	elapsed = cputime-tstart;
	d.timings{iN}{iR}{evalId} = elapsed;
	% store result of last computation for comparison
	res1s{evalId} = res1; res2s{evalId} = res2; 
      end % for evalId
      
      % set not computed times to NaN
      if (N>200)
	d.timings{iN}{iR}{1} = NaN;
	res1s{1} = res1s{3}; res2s{1} = res2s{3}; 
      end
      if (N>500)
	d.timings{iN}{iR}{2} = NaN;
	res1s{2} = res1s{3}; res2s{2} = res2s{3}; 
      end
      
      % compute the difference between the results
      errall = [
	  norm(res1s{1}-res1s{2},inf)
	  norm(res2s{1}-res2s{2},inf)
	  norm(res1s{1}-res1s{3},inf)
	  norm(res2s{1}-res2s{3},inf)
	  norm(res1s{2}-res1s{3},inf)
	  norm(res2s{2}-res2s{3},inf)
	       ];
      errtotal = max(errall);
      disp([errtotal errall']);
      %disp(['errtotal = ' num2str(errtotal)]);
      if (errtotal>1e-12)
	disp(errall)
	warning('difference in evaluation schemes larger than 1e-12.')
      end
    end % for iR
  end % for iN
  
  save('TwoD_CompareEvaluationSchemes.mat','d');
  rmpath('..')
  return
end


if strcmp(runOrPlot, 'plot')

  load('TwoD_CompareEvaluationSchemes.mat')

  set(0,'DefaultAxesFontSize',15)
  set(0,'DefaultTextFontSize',15)
  set(0,'DefaultLineLineWidth',1)
  set(0,'DefaultAxesLineWidth',1)
  set(0,'DefaultTextInterpreter','latex');
  figure(1)

  timesum = zeros(length(d.timings),1);
  timefft1 = timesum;
  timefft2 = timesum;
  for iR = 1:length(d.Rvals);
    for iN = 1:length(d.timings)
      timesum(iN) = d.timings{iN}{iR}{1};
      timefft1(iN) = d.timings{iN}{iR}{2};
      timefft2(iN) = d.timings{iN}{iR}{3};
    end
    %[timesum timefft1 timefft2]
    lh = loglog(d.Nvals,timesum, 'g-s', ...
		d.Nvals, timefft1, 'r-x', ...
		d.Nvals, timefft2, 'b-d');
    hold on
    %pause
  end
  c0 = 2*d.Rvals(2); 
  loglog(d.Nvals, 3e-6*d.c1*d.c2*c0^2*(d.Nvals.^4),'--', ...
	 'color', 0.7*[0.5 1 0.5],'linewidth',2)
  loglog(d.Nvals, 8e-7*d.c1^2*d.c2*d.Nvals.^3,'--', ...
	 'color', 0.7*[1 0.5 0.5],'linewidth',2);
  loglog(d.Nvals, 1e-7*d.c1*d.c2*4*(d.Nvals.^2).*log2(d.Nvals),'--', ...
	 'color', 0.7*[0.5 0.5 1],'linewidth',2);
  axis([4, 1000, 7e-3 1e4])
  xlabel('$N$'), 
  ylabel('CPU time (seconds)')
  
  [LEGH,OBJH,OUTH,OUTM] = ...
      legend(lh, {'summation', 'FFT1', 'FFT2'}, 'location', 'northwest');
  for hhi = 1:length(OBJH)
    hh = OBJH(hhi);
    if strcmp(get(hh,'type'), 'text')
      set(hh,'interpreter','latex');
    end
  end

  rmpath('..')
  return
end
