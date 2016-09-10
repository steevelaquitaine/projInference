

%slCalculateDiffAICmodelsvsSwitch.m
%
%
% author : steeve laquitaine
%purpose : calculate difference in AICs between a reference models ('switching')
%           and other models (Bayesian models)
%

function [AICsvsSw,semAICsvsSw,aicDiff] = slCalculateDiffAICmodelsvsSwitch(aics_ref,aics)

nModels =  size(aics,2);

aics_ref =  SLmakeColumn(aics_ref);

%average difference between models AICs and switch AIC over model
aic_ref_all = repmat(aics_ref,1,nModels);
aicDiff = aics - aic_ref_all;
AICsvsSw = nanmean(aicDiff,1);
semAICsvsSw = sem(aicDiff,1);
figure('color','w')
SLdrawBar(AICsvsSw,1:nModels,1:nModels,'yError',semAICsvsSw,'facecolor',[.5 .5 .5])
