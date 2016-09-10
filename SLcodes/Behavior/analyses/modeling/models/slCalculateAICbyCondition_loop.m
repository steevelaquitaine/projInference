


%slCalculateAICbyCondition_loop.m
%
%
% author: steeve laquitaine
%purpose: matrix of AICs by subject by model
%
%
%
%usage:
%         %by conditions for Bayes with unimodel peaks
%         models = {'model_CompDiv','model_Bayes_Sampling','model_Bayes_MAP'};
%         AICs = slCalculateAICbyCondition_loop(models);
%         [AICsvsSw00,semAICsvsSw00,aicDiff00] = slCalculateDiffAICmodelsvsSwitch(AICs(:,1),AICs(:,2:end))
%
%
%         %all data for others
%         model_ref = 'model_CompDiv';
%         models = {'model_Bayes_Sampling_withCard', ...
%             'model_Bayes_WJMtailedPrior','model_Bayes_WJM','model_Bayes_MAP_withCard',...
%             'model_Bayes_MAP_FatTailPrior'}; 
%         path = '/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/';
%         [AICsvsSw01,semAICsvsSw01,aicDiff01] = slPlotModelsAICvsSwitchingAIC(model_ref,models,path);
%
%         figure('color','w')
%         SLdrawBar([AICsvsSw00 AICsvsSw01],1:7,1:7,'yError',[semAICsvsSw00 semAICsvsSw01],'facecolor',[.5 .5 .5])

function AICs = slCalculateAICbyCondition_loop(models)

%loop over models
for i = 1 : length(models)
    path = ['/Volumes/DroboBKUP/data/dataPsychophy/proj01_priorStrength/modelfit/AIC/' models{i}];
    [~,allVars,sub] = slLoadSavedFitParams('allVars',path);
    AICs(:,i) = slCalculateAICbyCondition(allVars);
end
