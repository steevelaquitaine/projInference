
%slCalculateAICbyCondition.m
%
%
% author: steeve laquitaine
%purpose: calculate AIC for one model and selected sets of conditions
%
%usage:
% 
%        path = '/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/model_Bayes_MAP/';
%        [~,allVars,sub] = slLoadSavedFitParams('allVars',path)
%        AICs = slCalculateAICbyCondition(allVars)
% 
%

function AICs = slCalculateAICbyCondition(allVars)

%loop over subjects
for i = 1 : length(allVars)
    
    %select trials when motion directions was displayed 90 degrees away
    %from prior mean
    dist2Prior = SLvectors2signedAngle(allVars(i).disp,225,'polar');
    sortedTrials = abs(dist2Prior)>=90;    
    fprintf('%s \n','(slCalculateAICbyCondition) Trials with motion direction 90 deg away from priors have been selected')
                
    %calculate logl
    logl = sum(allVars(i).loglpertrial(sortedTrials));
    AICs(i) = 2*(allVars(i).nfitp - logl);    
end







