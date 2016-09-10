
%slfitLearningBayes.m
%
%
%
%    author: steeve laquitaine
%      date: 160608
%   purpose: fit learning Bayesian model with maximum likelihood fit with all parameters free
%
%     usage:
%
%        [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp,fithistory] = slfitLearningBayes(initp,data,disp,StimStrength,...
%            pstd,priorShape,priorModes,TheModel,options,varargin)
%
%nested functions : slfitLearningBayesSingleParamSet

function [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp,fithistory] = slfitLearningBayes(initp,data,disp,StimStrength,...
    pstd,priorModes,trials,TheModel,options,varargin)

%loop over sets of initial parameters
fprintf('%s \n','=================== fitting info ==============================')
fprintf('%s \n','(slfitLearningBayes) Algorithm: fminsearch-Nelder Mead')
fprintf('%s \n',['(slfitLearningBayes) init parameter set: ',...
    num2str(1),'/',num2str(size(initp,1))])
fprintf('%s \n','(slfitLearningBayes) priors are fixed')
fprintf('%s \n','=================================================================')

%fit
[fitPtmp,negLogl,exitflag,outputFit,history] = slfitLearningBayesSingleParamSet(initp,data,disp,StimStrength,...
    pstd,priorModes,trials,TheModel,options,varargin{:});

%Fit parameters and SSE
fitPbkp(1,:) = fitPtmp;
negLoglbkp(1) = negLogl;
exitflagbkp{1} =  exitflag;
outputFitbkp{1} =  outputFit;
fithistory{1}.hist = history;
fithistory{1}.initp = initp;

output.fitAlgo = 'fminsearch';
save('fithistory','fithistory')

end

