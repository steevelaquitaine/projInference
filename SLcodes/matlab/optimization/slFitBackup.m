
%slFitBackup.m
%
% author: steeve laquitaine
%   date: 150822
%purpose: backup fminsearch fit iterations, objective function and parameters and save them in a fitBackup.mat file.
%
%  usage: insert in optimset
%         
%         function [history] = myfit
%           
%           %init stored path
%           mkdir(['fitBackup01])
%           cd(['fitBackup01])
%           
%           %init stored variables
%           history.params    = [];
%           history.iteration = [];
%           history.funccount = [];
%           history.fval      = [];
%           
%           options = optimset('OutputFcn', @slFitBackup);
%           [x fval] = fminsearch(@objfun, x0,options);
%
%          end
%
function stop = slFitBackup(x,optimvalues,state);

  stop = false;
  if isequal(state,'iter')
    history.params    = [history.params    ; x                    ];
    history.iteration = [history.iteration ; optimvalues.iteration];
    history.funccount = [history.funccount ; optimvalues.funccount];
    history.fval      = [history.fval      ; optimvalues.fval     ];
  end

  save('fitBackup','history')
end