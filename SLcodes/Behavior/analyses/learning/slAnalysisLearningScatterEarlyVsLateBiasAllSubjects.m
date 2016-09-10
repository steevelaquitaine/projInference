

%slAnalysisLearningScatterEarlyVsLateBiasAllSubjects.m
%
%
% author: steeve laquitaine
%purpose: Scatter plot the slopes of estimates vs true motion direction for
%         early versus late trials pooled across runs. 
%         see "slAnalysisLearningScatterEarlyVsLateBias.m"
%
%  usage:
%   
%   [linfit_early,linfit_late,r,p] = slAnalysisLearningScatterEarlyVsLateBiasAllSubjects({'sub01',...
%       'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09',...
%       'sub10','sub11','sub12'},'~/data/dataPsychophy/proj01_priorStrength/','vonMisesPrior')
%
%
%Inputs:
%
%    subjects : 
%         e.g., {'sub01','sub02'}
%
%    dataPath :
%         e.g., '~/data/dataPsychophy/proj01_priorStrength/'
%     
%    'experiment': 
%          e.g. 'experiment','vonMisesPrior'
%
%Outputs:
%
%


function [linfit_early,linfit_late,r,p,d] = slAnalysisLearningScatterEarlyVsLateBiasAllSubjects(subjects,dataPath,experiment)

%loop over subjects and scatter plots slopes 
%for early vs late
nSub = length(subjects);
%organize fig axes
nax = round(sqrt(nSub+4));
for i = 1 : nSub
  
  %show all subject in a single figure
  newfig = figure(10000);
  set(gcf,'color','w')
  newAx = subplot(nax,nax,i); 
    
  %database
  d{i} = SLMakedatabank(subjects(i),'dataPath',dataPath,'experiment',experiment);
  
  %scatter plot
  [linfit_early{i},linfit_late{i},r{i},p{i}] = slAnalysisLearningScatterEarlyVsLateBias(d{i},100);
  
  %close variability figure, not informative
  close(figure(2))
  title([subjects{i} ': correlation: R=' num2str(round(r{i}(2,1)*100)/100) '; p=' num2str(p{i}(2,1))])
  set(gcf,'outerposition',[2540 752 411 393])
  
  %copy this plot to all subject plot  
  mygca = gca;
  slgraphCopyFigAxtoNewfigAx(mygca,newfig,newAx)    
end
