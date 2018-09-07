function fname = OneD_ConvergenceNonLocalTermApprox(R, OmegaId, sink, cosk, NR)


if (length(sink)>0)
  sinstr = ['_sin' strrep(num2str(sort(sink),'%2d'),' ','sin')];
else
  sinstr = '';
end
if (length(cosk)>0)
  cosstr = ['_cos' strrep(num2str(sort(cosk),'%2d'),' ','cos')];
else
  cosstr = '';
end
fname = ['OneD_ConvergenceNonLocalTermApprox' ...
	 '_R' strrep(num2str(R),'.','p') ... 
	 '_Omega' num2str(OmegaId) ...
	 sinstr cosstr ...
	 '_NR' num2str(NR,'%04d')
	];

oldpath = path;
addpath('..');


set(0,'DefaultAxesFontSize',20)
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextInterpreter','latex');

linetypes1 = {'r-x','g-s', 'b-d','c-+', 'y-o','k-v'};
linetypes2C = {'b.', 'g.', 'r.', 'c.', 'm.', 'y.', 'k.'};
linetypes2BW = {'k-x', 'k:', 'k-.', 'k--', 'k-'};

linetypes2 = linetypes2BW;


evalId  = 2; % 2=FFT
%Nvals = [50 100 200:200:2000 3000:1000:2000]+5;
Nvals = round(1./logspace(log10(1/55),log10(1/2005),20))
plotiN = [1:4:length(Nvals)]; % plot only every fourth value in
                              % spatially resolved error plot 
ruleIdvals=[1,2,3,4,6];

legstr={};
for iN = 1:length(Nvals)
  % setup spatial grid
  N = Nvals(iN);
  disp(['Using N = ' num2str(N)]);
  h = 1/N;
  N1=2*N;
  x=([0:(N1-1)]'+0.5)*h;
  
  % evaluate integral analytically
  g = zeros(size(x));
  aexact = g;
  
  for kk=1:length(sink)
     g = g + sin(2*pi*sink(kk)*x);
     aexact = aexact + sol1D_Fourierg(x+0.5*h, R, sink(kk), true, OmegaId);
  end
  for kk=1:length(cosk)
     g = g + cos(2*pi*cosk(kk)*x);
     aexact = aexact + sol1D_Fourierg(x+0.5*h, R, cosk(kk), false, OmegaId);
  end
  % constant g
  %g = ones(size(g));
  %aexact = 0*g;
  % linear g
  %g = x*7;
  %aexact = ones(size(g))*1/2*7; % in the interior of the domain for Omega_1 
  
  
  % evaluate non-local term using various methods to compute the weights
  for iruleId = 1:length(ruleIdvals)
    ruleId = ruleIdvals(iruleId);
    mask = setupIntegralRule1D(R, h, N1, OmegaId, ruleId, NR);
    a{ruleId} = evalIntegral1D(g, mask, evalId);
    % compute infty-norm of error
    errdiff=a{ruleId}-aexact;
    errinf{ruleId}(iN)=norm(errdiff, inf);
    errl2{ruleId}(iN) =norm(errdiff, 2)/sqrt(length(errdiff));
    % plot absolute errors 
    if any(iN==plotiN)
      disp(['Plotting in figure ' num2str(ruleId)])
      if (iN==1)
	legstr{iN}=['$N=' num2str(N) '$'];
      else
	legstr{iN}=['$' num2str(N) '$'];
      end
      figure(ruleId);
      legh{ruleId}(find(iN==plotiN)) = ...
	  semilogy(x, abs(errdiff), linetypes2{find(iN==plotiN)});
      hold on
      axis([x(1)-h/2 x(end)+h/2 1e-7 1e-1])
    end
    if (iN == plotiN(end))
      figure(ruleId);
      xlabel('$x$'),ylabel('absolute error')
      %legend(legh{ruleId}(plotiN),legstr{plotiN},'orientation', ...
      %       'horizontal');
      [LEGH,OBJH,OUTH,OUTM] = legend(legh{ruleId}(:),legstr{plotiN},...
				     'orientation', 'horizontal'); 
      for hhi = 1:length(OBJH)
	hh = OBJH(hhi);
	if strcmp(get(hh,'type'), 'text')
	  set(hh,'interpreter','latex');%,'FontSize',15);
	end
      end
    end
    drawnow;
  end % for iruleId = 1:length(ruleIdvals)
end % for iN = 1:length(Nvals)

% ...and plot spatially not resolved
format short e
[1./Nvals' ...
 errinf{1}(:) errinf{2}(:) errinf{3}(:) errinf{4}(:) errinf{6}(:)
]
figure(10)
for iruleId = 1:length(ruleIdvals)
  ruleId = ruleIdvals(iruleId);
  loglog(1./Nvals, errinf{ruleId}(:), linetypes1{ruleId});
  hold on
end
axv = axis;
axis([1/Nvals(end)*(1-1/10) 1/Nvals(1)*(1+1/10) 1e-5 1e-1])
% legend
legstr = {'(E-c)', '(E-l)', '(M-c)', '(M-l)', '(T-l)'};
[LEGH,OBJH,OUTH,OUTM] = legend(legstr, 'location','SouthEast'); 
for hhi = 1:length(OBJH)
  hh = OBJH(hhi);
  if strcmp(get(hh,'type'), 'text')
    set(hh,'interpreter','latex');
  end
end

%order slope
loglog(1./Nvals, 1500*(1./Nvals).^2,'--', ...
       'color',0.5*[1 1 1],'linewidth', 2);
xlabel('$h$')
ylabel('$\left\|\underline{a}\textrm{--}\mathcal{A}\{g(\cdot)\}\right\|$')

	 
rmpath('..')
return
