function TwoD_ConvergenceNonLocalTermApprox(runOrPlot)

%{
  diary('TwoD_ConvergenceNonLocalTermApprox.diary');
  TwoD_ConvergenceNonLocalTermApprox('run');
  TwoD_ConvergenceNonLocalTermApprox('plot')
  printeps('TwoD_ConvergenceNonLocalTermApprox.eps');
  diary off
%}


format compact
format short e


oldpath = path;
addpath('..');

if strcmp(runOrPlot, 'run')
  R=0.1;
  OmegaId = 2; % 1 or 2
  ruleId  = 1; % 1 or 2
  evalId = 3;  % FFT2

  k1vals = [1 2 4 8];
  k2vals = [1 2 1 1];
  Nvals = [round(1./logspace(log10(1/10),log10(1/500),10)) 1000]

  % store data in results structure
  res.R       = R;
  res.OmegaId = OmegaId;
  res.ruleId  = ruleId;
  res.evalId  = evalId;
  res.k1vals  = k1vals;
  res.k2vals  = k2vals;
  
    
  for iN = 1:length(Nvals)
    N = Nvals(iN);
    disp(['(N,R) = ' num2str([N,R])]);
    
    % setup spatial grid for function g
    h = 1/N;
    N1=N;
    N2=N;
    x=([0:(N1-1)]'+0.5)*h;
    y=x;
        
    % compute right-hand sides g
    [xx,yy] = meshgrid(x,y);
    gmat = zeros(prod(size(xx)),length(k1vals));
    for ik=1:length(k1vals)
      g = sin(2*k1vals(ik)*pi*xx).*cos(2*k2vals(ik)*pi*yy);
      gmat(:,ik) = g(:);
    end

    % select Omega function
    switch OmegaId
     case 1
      nl.OmegaStr = '@(r)(1/(pi*R^2))';
     case 2
      nl.OmegaStr = '@(r)(3/(pi*R^2)*(1-r/R))'; 
     otherwise
      error('unknown case.')
    end
    
    % compute integration weights
    % set initial value for number of subdivision
    if (iN ==1) % set initial value for number of subdivision
      sd = 20; 
    else
      sd=max(10, round(sd/2));
    end
    acc = max(min(h^2, 1e-3), 1e-5); % set desired accuracy of weights
    nl.x1Dir = setupIntegralRule2D('x1', R, h, N1, N2, eval(nl.OmegaStr), ...
					 ruleId, sd, acc);
    sd = nl.x1Dir.NR;
    nl.x2Dir = setupIntegralRule2D('x2', R, h, N1, N2, eval(nl.OmegaStr), ...
					 ruleId, sd, inf);
    disp('... rule setup complete')
    
    % evaluate nonlocal term
    a1mat = zeros(size(gmat));
    a2mat = zeros(size(gmat));
    for ik=1:length(k1vals)
      [a1, a2] = evalIntegral2D(gmat(:,ik), nl, evalId);
      a1mat(:,ik) = a1(:);
      a2mat(:,ik) = a2(:);
    end
    
    % save results
    res.N = N;
    res.h = h;
    res.a1mat = a1mat;
    res.a2mat = a2mat;
    res.x = x;
    res.y = y;
    
    fname = ['TwoD_ConvergenceNonLocalTermApprox_OmegaId' num2str(OmegaId) ...
	     '_ruleId' num2str(ruleId) ...
	     '_N' num2str(N,'%05d') '.mat'];
    disp(['Saving to file ' fname]);
    save(fname,'res');
    
    res.a1mat = [];
    res.a2mat = [];
    res.x = [];
    res.y = [];
    res.N = [];
    res.h = [];

  end %for iN = 1:length(Nvals)

  rmpath('..')
  return
end

if strcmp(runOrPlot, 'plot')
  
  OmegaId = 2; % 1 or 2
  ruleId  = 1; % 1 or 2
  Nvals = [round(1./logspace(log10(1/10),log10(1/500),10)) 1000]
  Nref = Nvals(end);
  Nvals = setdiff(Nvals,Nref); % remove Nref from Nvals if in there
  
  set(0,'DefaultAxesFontSize',15)
  set(0,'DefaultLineLineWidth',1)
  set(0,'DefaultAxesLineWidth',1)
  set(0,'DefaultTextFontSize',15)
  set(0,'DefaultTextInterpreter','latex');
  figure
  
  
  % load reference solution
  fname = ['TwoD_ConvergenceNonLocalTermApprox_OmegaId' num2str(OmegaId) ...
	   '_ruleId' num2str(ruleId) ...
	   '_N' num2str(Nref,'%05d') '.mat'];
  disp(['Loading reference solution from ' fname]);
  load(fname,'res');
  ref = res;
  res = [];
  [XX,YY] = meshgrid(ref.x, ref.y);
  for ik = 1:length(ref.k1vals)
    a1ref{ik} = reshape(ref.a1mat(:,ik), length(ref.y), length(ref.x));
    a2ref{ik} = reshape(ref.a2mat(:,ik), length(ref.y), length(ref.x));
  end
  
  % compute errors
  errmat=zeros(length(ref.k1vals),length(Nvals));
  
  legstr = {};
  for iN = 1:length(Nvals)
    N=Nvals(iN);
    % load solution
    fname = ['TwoD_ConvergenceNonLocalTermApprox_OmegaId' num2str(OmegaId) ...
	     '_ruleId' num2str(ruleId) ...
	     '_N' num2str(N,'%05d') '.mat'];
    disp(['Loading solution from ' fname]);
    load(fname,'res');
 
    % compute values to compare with from reference solution
    [xx,yy] = meshgrid(res.x, res.y);
    cmpa1=zeros(size(res.a1mat));
    cmpa2=cmpa1;
    for ik = 1:length(res.k1vals)
      % make sure here that the x and y coordinates are correct
      % (must be on the center of cell faces because we compute the
      % non-local term there!)
      a1 = interp2(XX+ref.h/2,YY,a1ref{ik},xx+res.h/2,yy,'cubic');
      cmpa1(:,ik) = a1(:);
      a2 = interp2(XX,YY+ref.h/2,a2ref{ik},xx,yy+res.h/2,'cubic');
      cmpa2(:,ik) = a2(:);
      legstr{ik}=['$(k_1,k_2)=(' num2str(res.k1vals(ik)) ',' ...
		  num2str(res.k2vals(ik)) ')$'];
    end
    
    % compute maximum absolute difference for each of the g-functions
    diff1 = cmpa1-res.a1mat;
    diff2 = cmpa2-res.a2mat;
    for ik = 1:length(res.k1vals)
      errmat(ik,iN) = max(norm(diff1(:,ik),inf), norm(diff2(:,ik),inf));
    end
  end
 
  [Nvals; errmat]
  linetypes1 = {'r-x','g-s', 'b-d','c-+', 'y-o','k-v'};
  for ik = 1:length(res.k1vals)
    lh(ik) = loglog(1./Nvals, errmat(ik,:),linetypes1{ik});
    hold on
  end
%  loglog(1./Nvals, errmat','-x')
  [LEGH,OBJH,OUTH,OUTM] = legend(lh,legstr,'location','northwest');
  for hhi = 1:length(OBJH)
    hh = OBJH(hhi);
    if strcmp(get(hh,'type'), 'text')
      set(hh,'interpreter','latex');%,'FontSize',15);
    end
  end
  loglog(1./Nvals, 1e3*(1./Nvals).^2,'--', ...
         'color',0.5*[1 1 1],'linewidth', 2); % second order slope
  xlabel('$h$')
  ylabel('$\left\|\underline{a}\textrm{--}\mathcal{A}\{g(\cdot)\}\right\|$')
  axis([1e-3 1.5e-1 1e-5 2])

  
  rmpath('..')
  return
end


