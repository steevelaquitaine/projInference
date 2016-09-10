
%slfitBayesfminSlogl.m
%
%
%
%    author: steeve laquitaine
%   purpose: Maximum likelihood fit with all parameters free


function [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp] = slfitBayesfminSlogl(initp,data,disp,StimStrength,...
    pstd,priorShape,priorModes,TheModel,options,varargin)


%loop over sets of initial parameters
parfor i = 1 : size(initp,1)
    
    %status
    fprintf('%s \n','=================== fitting info ==============================')
    fprintf('%s \n','(slfitBayesfminSloglPriorsFixed) Algorithm: fminsearch-Nelder Mead')
    fprintf('%s \n',['(slfitBayesfminSloglPriorsFixed) init parameter set: ',...
        num2str(i),'/',num2str(size(initp,1))])
    fprintf('%s \n','(slfitBayesfminSloglPriorsFixed) priors are fixed')
    fprintf('%s \n','=================================================================')
    
    %Nelder-Mead
    [fitPtmp,negLogl,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
        SLgetLoglBayesianModel(data,...
        disp,...
        StimStrength,...
        pstd,...
        fitPtmp,...
        priorShape,...
        priorModes,...
        TheModel,varargin{:}),...
        initp(i,:),...
        options);
    
    %Fit parameters and SSE
    fitPbkp(i,:) = fitPtmp;
    negLoglbkp(i) = negLogl;
    exitflagbkp{i} =  exitflag;
    outputFitbkp{i} =  outputFit;
end
output.fitAlgo = 'fminsearch';