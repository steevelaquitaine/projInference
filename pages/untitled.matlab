

%Basic Bayes 
%-----------
%datapath = '~/Desktop/dataPsychophy/proj02_priorStrength80';
%output = SLfitBayesianModel({'sub04'},[12 3 1 0.5 NaN NaN NaN NaN 0.03 15 NaN],'pscaling', ones(1,11),'dataPathVM',datapath,'experiment','vonMisesPrior','filename','datafit','MAPReadout','MaxLikelihoodFit', 'fminsearch'); 

%Switching like or prior
%-----------------------
% datapath='~/Desktop/dataPsychophy/proj02_priorStrength80'
% output = SLfitCompetitionModel({'sub04'},[25 5 1.4 0.2 NaN NaN NaN NaN 0.1 15], 'pscaling',ones(1,10),'dataPathVM',datapath,'experiment','vonMisesPrior', 'MaxLikelihoodFit','fminsearch', 'filename','datafit04sw','savedfolder','SwitchingLikePrior');

%Swicth post or prior
%--------------------
%datapath='~/Desktop/dataPsychophy/proj02_priorStrength80'
%output = SLfitCompetitionModel({'sub04'},[15 4 1 0.2 NaN NaN NaN NaN 0.1 30], 'pscaling',ones(1,10), 'dataPathVM',datapath,'experiment','vonMisesPrior', 'MaxLikelihoodFit','fminsearch', 'filename','datafitswPoOrPr80sub04', 'switchPostorPrior','savedfolder','SwitchingPriorPost')

%Basic Bayes tailed prior 
%------------------------
% datapath = '~/Desktop/dataPsychophy/proj02_priorStrength80'
% output = slshortfitBayesTailedPrior({'sub04'},[12 3 1 0.5 NaN NaN NaN NaN 0.03 15 0],ones(1,11),datapath);

%Sampling
%---------
% addpath(genpath('~/proj/steeve/SLcodes'))
% addpath(genpath('~/proj/mgl'))
% addpath(genpath('~/proj/mrTools'))
% datapath ='/Volumes/MoonData/data/dataPsychophy/proj02_priorStrength80'
% output = slshortfitBayesSmp({'sub04'},[20 5 2 0.8 NaN NaN NaN NaN 0.07 60 NaN],ones(1,11),datapath);

%Bayes cardinal
%--------------
% addpath(genpath('~/proj/steeve/SLcodes'))
% addpath(genpath('~/proj/mgl'))
% addpath(genpath('~/proj/mrTools'))
% datapath ='~/Desktop/dataPsychophy/proj02_priorStrength80'
% output = SLfitBayesianModel({'sub04'},[12 3 1 0.5 NaN NaN NaN 0 0.03 15 NaN],'pscaling',ones(1,11),'dataPathVM',datapath,'experiment','vonMisesPrior','filename','datafit','MAPReadout','MaxLikelihoodFit', 'fminsearch');

%Bayes sampling card
%-------------------
%datapath ='/Volumes/MoonData/data/dataPsychophy/proj02_priorStrength80'
%output = slshortfitBayesSmp({'sub04'},[20 5 2 0.8 NaN NaN NaN 0 0.07 60 NaN],ones(1,11),datapath);