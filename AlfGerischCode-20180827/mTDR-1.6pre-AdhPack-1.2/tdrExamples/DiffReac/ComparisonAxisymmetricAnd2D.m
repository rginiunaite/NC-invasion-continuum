function ComparisonAxisymmetricAnd2D(selectICandBC, reactiontype)
% function ComparisonAxisymmetricAnd2D(selectICandBC, reactiontype)
%
% This function is part of the DiffReac TDR model.
%
% This examples compares the effect on the solution of a reaction-diffusion
% equation caused by a standard 2D domain vs. an axisymmetric domain. 
%
%  selectICandBC ... integer, allowed are values 1 and 2
%  reactiontype  ... integer, allowed are values 1, ..., 6
%
% Examples:
%   ComparisonAxisymmetricAnd2D(1, 2);
%   ComparisonAxisymmetricAnd2D(1, 5);
%
% The full model consists of the six m-files: ComparisonAxisymmetricAnd2D.m,
% ProbBCs.m, ProbFReac.m, ProbFTrans.m, ProbFy0.m, and ProbGetParams.m .
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrExamples/DiffReac/ComparisonAxisymmetricAnd2D.m
%* Date created  : 2006, January 25
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     :
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2006 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************


clear functions
format compact;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% defaults for plots %%%
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultTextFontSize',15)
%%%%%%%%%%%%%%%%%%%%%%%%%%


oldpath=path();          % save current Matlab path
% include ./tdr and ./tdrUtil in path
path('../../tdrUtil',path); 
path('../../tdr',path); 


% print TDR directory and problem directory used 
% as convenience for the user
which tdrInit            % is in tdr directory
which ProbGetParams      % is in problem directory
    

% set param and options for ODE solver as empty
param = [];
options = [];
% install output function for odesolver
%options = odeset(options, 'OutputFcn', @outputfun);
% set tolerance options
options = odeset(options, 'AbsTol', 1e-8, 'RelTol', 1e-8);

% Select solution (valid are 1 or 2)
inputStruct.params.selectICandBC = selectICandBC;

% Select reaction term and set reaction term parameters
switch reactiontype
 case 1
  inputStruct.params.selectReaction = 'none'; 
  inputStruct.params.ReactionParams = []; % no params
  inputStruct.tdr.tvec = [0:100]*1/100;
 case 2
  inputStruct.params.selectReaction = 'logistic'; 
  inputStruct.params.ReactionParams = [1000 1];   % 2 params [alpha beta]
  inputStruct.tdr.tvec = [0:100]*0.02/100;
 case 3
  inputStruct.params.selectReaction = 'logistic'; 
  inputStruct.params.ReactionParams = [100 1];    % 2 params [alpha beta]
  inputStruct.tdr.tvec = [0:100]*0.1/100;
 case 4
  inputStruct.params.selectReaction = 'logistic'; 
  inputStruct.params.ReactionParams = [10 1];     % 2 params [alpha beta]
  inputStruct.tdr.tvec = [0:100]*0.5/100;
 case 5
  inputStruct.params.selectReaction = 'genlogistic'; 
  inputStruct.params.ReactionParams = [30 0.6];   % 2 params [alpha beta] 
  inputStruct.tdr.tvec = [0:100]*2/100;
 case 6
  inputStruct.params.selectReaction = 'genlogistic'; 
  inputStruct.params.ReactionParams = [100 0.6];  % 2 params [alpha beta]
  inputStruct.tdr.tvec = [0:100]*0.5/100;
 otherwise
  error('unknown reactiontype.')
end

% Select smallest y-coordinate of domain (>=0)
inputStruct.params.yshift = 0.01;

clear global TDRP
%
% RUN FOR STANDARD 2D DOMAIN
%
global TDRP;
timerReset();           % reset all timer
% init the TDR problem (initialises the data in TDRP and returns 
% the initial value y0 and the vector of output times tspan)
inputStruct.grd.isAxiSymmetric = false;
[y0, tspan] = tdrInit(inputStruct);
disp(['TDRP.grd.isAxiSymmetric = ' num2str(TDRP.grd.isAxiSymmetric)]);
% call to ODE solver performs the time integration
[t,y] = ode15s(@tdrFdgl, tspan, y0, options, param);
timerPrint('main::');   % print current values of timer
% Save one solution line parallel to the y-axis for all times
x2d   = TDRP.grd.y0(1)+ ([1:TDRP.grd.ny(1)]-0.5)* TDRP.grd.dy(1);
sol2d = y(:, 5*TDRP.grd.ny(1) + [1:TDRP.grd.ny(1)]);

clear global TDRP
%
% RUN FOR AXISYMMETRIC DOMAIN
%
global TDRP;
timerReset();           % reset all timer
% init the TDR problem (initialises the data in TDRP and returns 
% the initial value y0 and the vector of output times tspan)
inputStruct.grd.isAxiSymmetric = true;
[y0, tspan] = tdrInit(inputStruct);
disp(['TDRP.grd.isAxiSymmetric = ' num2str(TDRP.grd.isAxiSymmetric)]);
% call to ODE solver performs the time integration
[t,y] = ode15s(@tdrFdgl, tspan, y0, options, param);
timerPrint('main::');   % print current values of timer
% Save one solution line parallel to the y-axis for all times
xcyl  = TDRP.grd.y0(1)+ ([1:TDRP.grd.ny(1)]-0.5)* TDRP.grd.dy(1);
solcyl = y(:, 5*TDRP.grd.ny(1) + [1:TDRP.grd.ny(1)]);


%
% Make the comparison
%

if (selectICandBC ~= 1)
  error('Comparison only tested for selectICandBC == 1.');
end

figure;
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

ax(1)=axes('position',[0.09 0.537  0.86 0.39]);
ax(2)=axes('position',[0.09 0.297  0.86 0.2]);
ax(3)=axes('position',[0.09 0.057  0.86 0.2]);

axes(ax(1))
plot(t, solcyl(:,1), 'k-', t, sol2d(:,1), 'k--')
disp(['solutions, ' inputStruct.params.selectReaction ', '...
       num2str(inputStruct.params.ReactionParams)]);

t0 = inputStruct.tdr.tvec(1);
tend = inputStruct.tdr.tvec(end);
tspan = tend-t0;

axislimits=axis();
ypos=axislimits(4)*1.1;  % for windows matlab
ypos=axislimits(4)*1.07; % for linux matlab
xpos= t0+tspan/2;
switch reactiontype
 case {1}
  tstr = ['Reaction term none'];
 case {2,3,4}
  tstr = ['Logistic reaction term with $\alpha = ' ...
	  num2str(inputStruct.params.ReactionParams(1)) '$ and $\beta = ' ...
	  num2str(inputStruct.params.ReactionParams(2)) '$.'];
 case {5,6}
  tstr = ['Generalised logistic reaction term with $\alpha = ' ...
	  num2str(inputStruct.params.ReactionParams(1)) '$ and $\beta = ' ...
	  num2str(inputStruct.params.ReactionParams(2)) '$.'];
 otherwise
  error('unknown reactiontype.')
end  
text(xpos, ypos, tstr, 'HorizontalAlignment', 'center', ...
     'interpreter', 'latex');

%axis([t(1) t(end) 1e-2 1])
lh(1)=legend('a', 'b', 'Location', 'West');
set(lh(1),'Interpreter', 'latex', 'Box', 'off', 'FontSize', 24);
set(lh(1), 'String', {'$\tilde u(t,r_0)$', '$\bar u(t,r_0)$'});

axes(ax(2)); %%% title('absolute difference')
plot(t, -(solcyl(:,1)-sol2d(:,1)),'k-' )
%axis([t(1) t(end) 1e-5 1])
labelstring = ...
    '$\bar u(t,r_0) - \,\tilde u(t,r_0)$';
%text(0.0005, -0.02, labelstring, 'Interpreter', 'latex')

axes(ax(3)); %%% title('relative difference in percent')
plot(t, 100*(-(solcyl(:,1)'-sol2d(:,1)'))./max([abs(sol2d(:,1)'); ...
            abs(solcyl(:, 1)')]),'k-' )
labelstring = ['$\displaystyle\frac{\bar u(t,r_0) - ' ...
	       '\,\tilde u(t,r_0)}{\max\{|\bar u(t,r_0)|, |\tilde' ...
	       ' u(t,r_0)|\}}\cdot \!100$'];
%text(0.0005, -17, labelstring, 'Interpreter', 'latex')

for i=1:2
  set(ax(i),'XTickLabel', {})
end
for i=3:3
  set(ax(i),'XTickLabel', ...
	    {num2str(t0); num2str(t0+tspan/4); num2str(t0+tspan/2); ...
	     'time t'; num2str(tend)}, ...
	    'XTick', t0 + [0:4]*tspan/4);  
end

switch reactiontype
 case {1}
  filename = ['EPS/Reac_' inputStruct.params.selectReaction '.eps'];
 case {2,3,4,5,6}
  filename = ['EPS/Reac_' inputStruct.params.selectReaction '_' ...
	      num2str(inputStruct.params.ReactionParams(1)) '_' ...
	      num2str(inputStruct.params.ReactionParams(2)) '.eps'];  
 otherwise
  error('unknown reactiontype.')
end  
%disp(['saving to ' filename]);
%print('-deps2', '-loose', filename);


%
% reset Matlab path
path(oldpath);
  
return

% Old code to create a movie. Needs adjustment!!
if (0)
%
% PLOT COMPARISON
%
makemovie = false;
%makemovie = true;

f1=figure('PaperUnits','centimeters', 'PaperOrientation', 'landscape', ...
      'PaperType', 'A4');
set(f1,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
           'NextPlot','replace','Visible','off')
if (makemovie)
  movieFile = ['VergleichAxiUnd2D_ICandBC' num2str(inputStruct.params.selectICandBC) ...
           '.avi']; 
  mov = avifile(movieFile, 'FPS', 5, 'Compression', 'Indeo3', 'Quality', 90);
end

for ti=1:length(t)
  tiString = num2str(t(ti))
  subplot(121)
  plot(x2d, sol2d(ti,:),'r', xcyl, solcyl(ti,:), 'b');
  %legend('1D', 'radial');
  
  ymin = min([sol2d(ti,:) solcyl(ti,:)]');
  ymax = max([sol2d(ti,:) solcyl(ti,:)]');
  %axis([xcyl(1)-0.1 xcyl(end)+0.1 ymin*0.98 ymin + 1.3*(ymax-ymin)]);
  
  if (inputStruct.params.selectICandBC == 1)
    if (t(ti)>0.4)
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.5 1.0]);
    elseif (t(ti) > 0.2)
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.2 1.0]);
    else
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.0 1.0]);
    end
  elseif (inputStruct.params.selectICandBC == 2)
    if (t(ti)>0.06)
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.0 0.4]);
    elseif (t(ti) > 0.01)
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.0 0.6]);
    else
      axis([xcyl(1)-0.1 xcyl(end)+0.1 0.0 1.0]);
    end    
  else
    error('');
  end
  xlabel('Radial direction r');
  text('String',['\quad$\mathbf{\tilde u(t,r)}$\quad'], ...
       'Interpreter', 'latex', 'Color', 'b', 'Rotation', 0, ...
       'Units', 'normalized', 'Position', [0.5, 1.05], 'HorizontalAlignment', 'left')
  text('String',[' , '], ...
       'Interpreter', 'none', 'Color', 'k', 'Rotation', 0, ...
       'Units', 'normalized', 'Position', [0.5, 1.05], 'HorizontalAlignment', 'center',...
       'VerticalAlignment', 'middle')
  text('String',['\quad$\mathbf{v(t,r)}$\quad'], ...
       'Interpreter', 'latex', 'Color', 'r', 'Rotation', 0, ...
       'Units', 'normalized', 'Position', [0.5, 1.05], 'HorizontalAlignment', 'right')

  subplot(122)
  plot(x2d, 0*x2d, 'r');
  hold on
  % absolute difference
  errdata1 = sol2d(ti,:) - solcyl(ti,:);
  % relative difference
  errdata2 = errdata1./max([abs(sol2d(ti,:)); abs(solcyl(ti,:))])*100; 
  plot(x2d, errdata2, 'k');
  hold off

  if (inputStruct.params.selectICandBC == 1)
    axis([xcyl(1)-0.1 xcyl(end)+0.1 -1.5E-2 1E-2]);
    axis([xcyl(1)-0.1 xcyl(end)+0.1 -50 20]);
  elseif (inputStruct.params.selectICandBC == 2)
    axis([xcyl(1)-0.1 xcyl(end)+0.1 -1.5E-2 1E-2]);
    axis([xcyl(1)-0.1 xcyl(end)+0.1 -50 20]);
  else
    error('');
  end

  xlabel('Radial direction r');
  title(['$\mathbf{v(t,r) - \tilde u(t,r)}$'], ...
     'Interpreter', 'latex', 'Color', 'k')
  text('String', ['t = ' tiString], 'Position',  [0.3, 0.85], 'Units', 'normalized');
   
if (makemovie)
  F = getframe(gcf);
  mov = addframe(mov,F);
else
  pause %(0.01)
end

end

if (makemovie)
  mov = close(mov);
end

end % if (0,1)
