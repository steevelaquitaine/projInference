

%slGetMotionTaskCond.m
%
%author : steeve laquitaine
%usage:
%
%       [o,nCond] = slGetMotionTaskCond(stimFeatureDeg,StimStrength,pstd,priorModes,'vonMisesPrior')


function [o,nCond] = slGetMotionTaskCond(stimFeatureDeg,stimStrength,pstd,priorModes,priorShape)

%warning
if sum(strcmp(priorShape,'vonMisesPrior'))~=1
    if sum(strcmp(priorShape,'bimodalPrior'))~=1

        fprintf('(slGetMotionTaskCond) Please input a prior type \n')
        fprintf('(slGetMotionTaskCond) Your choices are:         \n')
        fprintf('(slGetMotionTaskCond) - "vonMisesPrior"         \n')
        fprintf('(slGetMotionTaskCond) - "bimodalPrior"          \n')
        dbstack
        keyboard
    
    end
    
%von Mises prior
elseif strcmp(priorShape,'vonMisesPrior')
    
    [o.uniqCond,~,o.posC] = SLuniqpair([pstd stimStrength stimFeatureDeg]);
    nCond = size(o.uniqCond,1);
    priorCond = pstd;
    numPriorCond = size(SLuniqpair(pstd),1); %check prior conditions

%bimodal prior 
elseif strcmp(priorShape,'bimodalPrior')
    
    priorModes = cell2mat(databank.priormodes); %prior
    priorCond = priorModes(:,2) - priorModes(:,1);
    numPriorCond = size(SLuniqpair(priorModes),1);%check prior conditions
end

if numel(unique(priorCond)) == numPriorCond
    
    [o.uniqCond,~,o.posC] = SLuniqpair([priorCond stimStrength stimFeatureDeg]);
    nCond = size(o.uniqCond,1);
    
end