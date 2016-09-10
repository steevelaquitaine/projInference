 %Author: Steeve Laquitaine
   %date: 140326
%purpose: lists analyses

%need to add the path of my code libary
addpath(genpath('/Users/steeve/Dropbox/Library/Library_code'))
%addpath(genpath('/Users/steeve_laquitaine/Dropbox/Library/Library_code'))


%% 140116 - draw HRF
drawHRF
set(gca,'xtick',cumsum([0.3 mean([1 2.3]) 3 0.1 0.1]),'xticklabel',...
    cumsum([0.3 mean([1 2.3]) 3 0.1 0.1]),...
    'fontsize',14)


%% 140117 - draw forward encoding model
figure('color','w');clf
simForwardEncoding

%% 140117 - make Bayesian inference predictions with forward encoding model
%figure('color','w')
clf
figure('color','w')
channel=simBayesWithForwardEncoding;

%% 140326 - Bayesian neural network (BNN)
%View firing rate and BOLD predictions for different motion directions
for motionDir=1:1:360;
    [MTBOLD,MTBOLDLIP,rLIPk,rMTj]=NNwithXGain(motionDir,'priorstd=0','numVoxelsMT=100','numVoxelsLIP=100',10,'display=on');
    drawnow;
end

%% 140326 - Simulate BOLD responses from a Bayesian neural network (training 
%and test data), run forward encoding model on training data and use the 
%channels weights to reconstruct motion directions in test data.
tic
SimulateFE00
toc


%% 140403 - Test predictions for areas encoding prior, llh and posterior 
%exclusively
tic
SimulateFE
toc
