function res = outputfunPlot(t, y, flag, params)
%
% function res = outputfunPlot(t, y, flag, params)
%
% This function can be used as output function of a Matlab ODE solver. 
% In this case it is installed with 
%   ODEoptions = odeset('OutputFcn', @outputfunPlot, [further options]);
% in struct ODEoptions, which is then supplied to the ODE solver.
% 
% This output function will then be called whenever the ODE solver 
% reaches one of the intermediate output time points. It opens figures
% one to tdr.size and plots the solution component of the TDR system 
% corresponding to the current time in each of these figures. The 
% minimum and maximum value of each component are displayed as well.
% Note that params is the variable param provided in the call to the ODE 
% solver and not the struct params generated in ProbGetParams().
%
%*********************  MATLAB TDR SYSTEM  ************************************
%* File          : tdrUtil/outputfunPlot.m
%* Date created  : 2005, September 29
%* Author(s)     : Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
%* Version       : 1.0
%* Revisions     :
%*
%*********************  COPYRIGHT NOTICE  *************************************
%* Copyright (C) 2004-2005 Alf Gerisch
%*                         Martin-Luther-University Halle-Wittenberg
%*                         Germany
%*
%* The TDR system in Matlab has been implemented by
%*   Mathias Franz (Oct 2004 - Feb 2005)
%*   Alf Gerisch   (Oct 2004 -         )
%******************************************************************************

global TDRP

if strcmp(flag, '')
  % output
  for EqNo = 1:TDRP.tdr.size
    figure(EqNo);
    clf;
    [ymin, ymax] =  tdrSolPlot(y(:,end), EqNo);
    title(['t=' num2str(t(end)) ',  EqNo = ' num2str(EqNo)]);
    xlabel('x'); ylabel('y'); zlabel('z');
    axis('equal');

    disp([t, EqNo, ymin,ymax])
    if (abs(ymin-ymax) < 1e-5)
      ymax  = ymin + 1e-5;
    end
    caxis([ymin,ymax]);

    colorbar
    hold('off')
    view(2) 
    drawnow
    %print('-depsc', ['output_t' num2str(t) '_Eq_' num2str(EqNo) '.eps']);
  end
end

res = 0;

return;
% end of function