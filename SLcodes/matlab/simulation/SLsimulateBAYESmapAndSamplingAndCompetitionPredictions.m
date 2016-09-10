
%   SLsimulateBAYESmapAndSamplingAndCompetitionPredictions.m
%
%     author: steeve laquitaine
%       date: 140601 updated 150109
%
%    purpose: simulate distribution pattern of motion direction estimates
%             predicted by Bayesian and Competition models (learnt priors,
%             cardinal priors, motor noise, random estimation).
%
%     usage:
%
%     SLsimulateBAYESmapAndSamplingAndCompetitionPredictions({'sub01'},...
%             [5 3 1 33 4.77 10.7 33 NaN 0 20],...
%            'dataPath','/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2');
% 
% 
%We test with reasonable parameters (similar to initial parameters obtained
%from manually fitting models to subjects data)
% - motor noise Km=20 (std=12 deg) 
% - fraction of andom estimation =0
%
%
%Typically the models predictions diverge most when the direction of the 
%displayed motion contradicts the prior and prior is relatively strong 
%(more competition)).


function SLsimulateBAYESmapAndSamplingAndCompetitionPredictions(subjects,simP,varargin)

%When 
%   - motion direction contradicts the prior
%   - evidence is weak and prior weak
%the bias/variability ratio is larger for MAP than for Sampling


%---------------------
%BAYES MAP predictions
%---------------------
%get already existing figures. We want to keep them.
figAlreadyOpened = SLgetFigures;
dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};

%simulate
SLsimulateBayesianModel(subjects,...
    simP,...
    'dataPath',dataPath,...
    'experiment','vonMisesPrior',...
    'MAPReadout');

%Copy figure we care about in a new figure
%weak llh and mild prior (1st axis) is the best condition
myAxes = SLgetAxes;
myAxes = myAxes(end:-1:1);
f1 = figure;
set(f1,'color','w')
ax1 = subplot(231);
copyobj(allchild(myAxes(1)),ax1);

%title
xlabel('MAPReadout')

%Keep only the figure we need and close unnecessary figures produced 
%by the previous simulation
figOpened = SLgetFigures;
figNewlyOpened = setdiff(figOpened,figAlreadyOpened);
close(figNewlyOpened(1:end-1))




%--------------------------
%BAYES Sampling predictions
%--------------------------
%check the figured opened outside the function
figAlreadyOpened = SLgetFigures;

%simulate
SLsimulateBayesianModel(subjects,...
    simP,...
    'dataPath',...
    dataPath,...
    'experiment','vonMisesPrior',...
    'SamplingReadout');

%Copy figure we care about in a new figure
%weak llh and mild prior (1st axis) is the best condition
myAxes = SLgetAxes;
myAxes = myAxes(end:-1:1);
figure(f1);
ax2 = subplot(232);
copyobj(allchild(myAxes(1)),ax2);

%title
xlabel('SamplingReadout')

%only keep the predictions figure
figOpened = SLgetFigures;
figNewlyOpened = setdiff(figOpened,figAlreadyOpened);
close(figNewlyOpened)




%-----------------------
%Competition predictions
%-----------------------
%check the figured opened outside the function
figAlreadyOpened = SLgetFigures;

%simulate
SLsimulateCompetitionModel(subjects,...
    simP,...
    'dataPath',...
    dataPath,...
    'experiment','vonMisesPrior');

%Copy figure we care about in a new figure
%weak llh and mild prior (1st axis) is the best condition
myAxes=SLgetAxes;
myAxes=myAxes(end:-1:1);
figure(f1);
ax3=subplot(233);   
copyobj(allchild(myAxes(1)),ax3);

%title
xlabel('Competition')

%only keep the predictions figure
figOpened=SLgetFigures;
figNewlyOpened=setdiff(figOpened,figAlreadyOpened);
close(figNewlyOpened)





