
%slPlotModelsAICvsSwitchingAIC
%
%
% author : steeve laquitaine
%purpose : load and plot difference in AICs between a reference models ('switching')
%           and other models (Bayesian models)
%
%   usage:
%
%         model_ref = 'model_CompDiv';
%         models = {'model_Bayes_Sampling','model_Bayes_Sampling_withCard', ...
%             'model_Bayes_WJMtailedPrior','model_Bayes_WJM','model_Bayes_MAP_withCard',...
%             'model_Bayes_MAP','model_Bayes_MAP_FatTailPrior'};
%         path = '/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/';
%         [AICsvsSw,semAICsvsSw,aicDiff] = slPlotModelsAICvsSwitchingAIC(model_ref,models,path)

%         model_ref = 'model_CompDiv';
%         models = {'model_Bayes_Sampling'};
%         path = '/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/';
%         [AICsvsSw,semAICsvsSw,aicDiff] = slPlotModelsAICvsSwitchingAIC(model_ref,models,path)


function [AICsvsSw,semAICsvsSw,aicDiff] = slPlotModelsAICvsSwitchingAIC(model_ref,models,path)

%load AICs
nModels = size(models,2);

%load model ref AIC
path_i = [path model_ref];
aics_ref = slLoadSavedFitParams('AIC',path_i);
    
%load other models AIC    
for i = 1 : nModels
    path_i = [path models{i}];    
    aics(:,i) = slLoadSavedFitParams('AIC',path_i);
end

%AIC differences
[AICsvsSw,semAICsvsSw,aicDiff] = slCalculateDiffAICmodelsvsSwitch(aics_ref,aics);




