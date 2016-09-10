

%slPlotPropSubSwitchingVsBayesModels.m
%
%
% author: steeve laquitaine
%purpose: load saved fit parameters for the different models of motion
%         direction estimation
%
%usage:
%
%        cd('~/proj/steeve/SLcodes/Behavior/data/)
%        [rowheader,colheader,data] = slcsvRead('modelfitData.csv');
%        [nRefWins,nNoWin,nRefLoses] = slPlotPropSubSwitchingVsBayesModels(rowheader,colheader,data,'Switching observer')

function [nRefWins,nNoWin,nRefLoses] = slPlotPropSubSwitchingVsBayesModels(rowheader,colheader,data,refmodel)

%get reference model data
refrow = find(strcmp(rowheader,refmodel));
refdata = data(refrow,:);

%get other models data
otherrows = setdiff(1:length(rowheader),refrow);
otherModelsData = data(otherrows,:);
othermodels = rowheader(otherrows,:);

nSubjects = length(colheader)-1;

for i = 1 : length(othermodels)    
     refLoses(i,:) = refdata - otherModelsData(i,:) > 2;
     noWin(i,:) = abs(refdata - otherModelsData(i,:)) < 2;
     refWins = 1 - (refLoses+noWin);
   
     %number of subjects
     nRefWins(i) = sum(refWins(i,:));
     nRefLoses(i) = sum(refLoses(i,:));
     nNoWin(i) = sum(noWin(i,:));
     
     %proportion of subjects for which reference model wins
     propRefWins(i,:) = sum(refWins(i,:))/nSubjects;
     propRefLoses(i,:) = sum(refLoses(i,:))/nSubjects;
     propNoWin(i,:) = sum(noWin(i,:))/nSubjects;
end

figure('color','w')
bar([propRefWins propNoWin propRefLoses]*100,'stacked')
colormap('gray')
box off
set(gca,'xtick',1:length(othermodels),'xticklabel',othermodels)
rotateXLabels(gca(),60)
ylabel('Reference model proportion of wins (%)')
title(['Reference model is: ' refmodel])
ylim([0 101])
