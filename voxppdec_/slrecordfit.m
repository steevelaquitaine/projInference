%slrecordfit.m
%
% author: steeve laquitaine
%purpose: record fit data in "history" variable
%
%
% usage:
%
%     global history
%     options.OutputFcn = slrecordfit;
%     [x objf] = fminsearch(@(x) exp(x(1))*(4*x(1)^2+2*x(2)^2+x(1)*x(2)+2*x(2)),[0 1],options)
% 
%     plot(history.iter,history.objf)

function outputfcn = slrecordfit

%Record fit values
global history
outputfcn = @record;
history.x    = [];
history.objf = [];
history.iter = [];

    function stop = record(x,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
            history.x = [history.x; x];
            history.objf = [history.objf; optimvalues.fval];
            history.iter = [history.iter; optimvalues.iteration];
        end
    end
end