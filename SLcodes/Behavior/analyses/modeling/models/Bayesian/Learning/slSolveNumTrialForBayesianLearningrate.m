

%slSolveNumTrialForBayesianLearningrate.m
%
%
% author: steeve laquitaine
%purpose: find trial number for 95% learning
%
%usage: 
%
%   t = slSolveNumTrialForBayesianLearningrate(0.55,0.95)
%
%note : the closed-form solution is t = -log(1-percLearnt)*tau

function t = slSolveNumTrialForBayesianLearningrate(tau,percLearnt)

%initialize trial number to 1
t0 = 1;

%re-write function so that it equals 0
percNotLearnt = percLearnt - 1;
fun = @(t) exp(-t/tau)+percNotLearnt;

%find trial number for 95% learning
t = fsolve(fun,t0);