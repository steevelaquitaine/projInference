
% SLsimulateCompetitionModel.m
%
%     author: steeve laquitaine
%       date: 140609 last modification 140611
%
%    purpose: simulate motion direction estimation behavior with
%             Competition between evidence and learnt prior.
%             Both representation can be biased by a built-in cardinal prior.
%             Sometimes subject makes random estimation.
%             Percepts are convolved with motor noise.
%
%      usage:
%%%lab
%
%unimodal
%--------
%sub09
% [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=SLsimulateCompetitionModel({'sub09'},...
%     [80 16 2 0.1 0.7 10 10 0 0.001 18],...
%     'dataPath',...
%     '/Users/steeve/Steeve_psychophy_data/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%     'experiment','vonMisesPrior');
%
%bimodal
%-------
% [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=SLsimulateCompetitionModel({'sub03'},...
%     [80 40 3 30 30 30 30 NaN 0.001 15],...
%     'dataPath',...
%     '/Users/steeve/Steeve_psychophy_data/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%     'experiment','bimodalPrior');
% 
%%mylaptop
%     [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=SLsimulateBayesianModel({'sub01'},...
%             [80 40 20 1.74 4.77 10.74 34.25 10.74 0.001 100],...
%            'dataPath',...
%            '/Users/steeve_laquitaine/Dropbox/Project_withJustin/data/psychophy_dataForPlayingAround/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%            'experiment','vonMisesPrior');
% 
%-------------
%varargins are
%--------------
%'dataPath',then 
%   path
%e.g., '/Users/steeve/Steeve_psychophy_data/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior',...
%
%
%'experiment',
%   - 'vomMisesPrior'
%   - 'bimodalPrior'
%
%
%description: we can simulate data with a Bayesian model 
%             (including learnt priors, a cardinal prior, 
%             motor noise and random choice)
%       
%             When we want to simulate fake data to check model fitting,
%             SlfitBayesianModel.m  can then be used to get back the 
%             model's parameters.
%
%             The model has 9 fit parameters
%               - 3 von Mises llh concentration parameters k (0<kl<inf, stimulus noise)
%               - 4 von Mises (or mixture of von Mises) priors concentration parameters 
%                 (0<kp<inf)
%               - 1 parameter for the fraction of trial with random estimation (0<prand<1)
%               - 1 parameter for motor noise (0<km<inf)
%
%             The model can use:
%               - von Mises learnt priors
%               - bimodal learnt priors (mixture of von Mises)
%
%             The model accounts for: 
%               - cardinal biases (Bayesian inference with cardinal priors)
%               - random estimation
%               - motor noise
%
%references: 
%     -Hurliman et al, 2002,VR
%     -Stocker&Simoncelli,2006,NN
%     -Girshick&Simoncelli,2011,NN
%     -Chalk&Series,2012,JoV
%
%
function [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=SLsimulateCompetitionModel(subjects,simP,varargin)

tic 

%call for help
if ieNotDefined('subjects')
    help SLsimulateCompetitionModel
  return    
end

%set data path
if sum(strcmp(varargin,'dataPath'))
    dataPath=varargin{find(strcmp(varargin,'dataPath'))+1};
    cd(dataPath)
else
    fprintf('%s \n','(SLsimulateCompetitionModel) You need to set dataPath')
   return
end

%Gather all available data in directory.
fprintf('%s \n','(SLsimulateCompetitionModel) Now creating a databank...')
databank=SLMakedatabank(subjects,varargin);

%Make predictions
%----------------
%   Data from each trial can be explained by either:
%
%       Bayesian inference with sensory evidence and cardinal
%       -----------------------------------------------------
%       measurement mi (=llh mean: ul) is drawn from a von Mises measurement 
%       density (v(mi,km)) in 1-("Pp"+"Prandom") fraction of trials.
%       is converted to likelihood and combined with learnt prior. The
%       probability of all possible percepts is the probability that it
%       comes from a given measuremen mi given by (v(mi,km)) * the
%       probability that subject switch given by a divisive competition
%       mechanism between evidence strength and learnt prior strength.
%       
%       Bayesian inference with learnt prior and cardinal prior
%       ---------------------------------------------------------
%       The probability of a percept given subject switched to the learnt
%       prior is the probability to switch given by the divisive
%       competition * the probability that this percept is produced by the
%       model(e.g., if cardinal is flat P(learnt prior mean)=100%;  or if 
%       learnt prior is flat P(each of 4 cardinal modes) = 25%)
%
%       Random estimation 
%       -----------------
%       an estimate is produced randomly in a fraction of trials ("Prandom").
%
%       And motor noise
%       ---------------
%       ~v(0,kmo) is added to the estimate (let's call it percept)
%       of those two processes.
%
%       Thus, the maximum likelihood of observing each trial data assuming this 
%       model is: 
%       probability of data given Bayesian inference * (1-"Pp"+"Prandom") + 
%       probability of data given random choice * "Prandom".
%       then both are convolved with motor noise (note: its is because convolution
%       is distributive).

%get experimental conditions
d=cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
coh=cell2mat(databank.data(:,(strcmp(databank.nm,'coh'))==1));
pstd=cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
priorModes=cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));

%subjects' data, %make sure 360 and 0 degrees are same,remove missing data
data=round(cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1)));
data(data==0)=360;
data=data(isnan(data)==0);
d=d(isnan(data)==0);
coh=coh(isnan(data)==0);
pstd=pstd(isnan(data)==0);

%get fit parameters. They should be in this order.
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal',
%'Prand','km'
kl1=simP(1);
kl2=simP(2);
kl3=simP(3);
klearnt1=simP(4);
klearnt2=simP(5);
klearnt3=simP(6);
klearnt4=simP(7);
kcardinal=simP(8);
Prandom=simP(9);
km=simP(10);
ulall=1:1:360;


%----------
%SET PRIORS
%----------
%
%(case Cardinal prior)
%---------------------
if isnan(kcardinal)
    %case we don't want to fit the cardinal prior.
    TheModel='withoutCardinal';
    kcardinal=0;
else
    %case we want to fit the cardinal prior.
    TheModel='withCardinal';
end

%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time.
%Constrain other fraction parameters between 0 and 1, and other parameters
%as>0.
%and case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom>1
    logL=inf;
    return
end
if any(simP<0)
    logL=inf;
    return
end
if any(isnan(simP([1,2,3,4,5,6,7,9,10])))
    logL=inf;
    fprintf('%s \n',['(getLogl) One of your fit parameter',... 
        ' that is not Kcardinal is NaN'])
    keyboard
end

%get prior mode(s)
SetOfpriorModes=SLuniqpair(priorModes);

%(case learnt prior is a von Mises)
%---------------------------------
if sum(strcmp(varargin,'vonMisesPrior'))==1
    
    priorShape='vonMisesPrior';
    
    %get prior mode (225 deg)
    modesPrior1=SetOfpriorModes(1,:);
    modesPrior2=SetOfpriorModes(1,:);
    modesPrior3=SetOfpriorModes(1,:);
    modesPrior4=SetOfpriorModes(1,:);
    
    %find trials of each prior
    Prior1=pstd==80;
    Prior2=pstd==40;
    Prior3=pstd==20;
    Prior4=pstd==10;
    
    %strength prior
    klearnt(pstd==80)=simP(4);
    klearnt(pstd==40)=simP(5);
    klearnt(pstd==20)=simP(6);
    klearnt(pstd==10)=simP(7);
end

%(case bimodal prior)
%--------------------
if sum(strcmp(varargin,'bimodalPrior'))
    priorShape='bimodalPrior';
    
    %get priors modes
    %[145 305]
    %165 285]
    %[185 265]
    %[205 245]
    modesPrior1=SetOfpriorModes(1,:);
    modesPrior2=SetOfpriorModes(2,:);
    modesPrior3=SetOfpriorModes(3,:);
    modesPrior4=SetOfpriorModes(4,:);

    %find the trials of each prior
    Prior1=priorModes(:,1)==modesPrior1(1,1) & priorModes(:,2)==modesPrior1(1,2);
    Prior2=priorModes(:,1)==modesPrior2(1,1) & priorModes(:,2)==modesPrior2(1,2);
    Prior3=priorModes(:,1)==modesPrior3(1,1) & priorModes(:,2)==modesPrior3(1,2);
    Prior4=priorModes(:,1)==modesPrior4(1,1) & priorModes(:,2)==modesPrior4(1,2);

    %strength (should be 20 deg)
    klearnt(pstd==20)=simP(6);
end

%-------
%set LLH
%-------
%mean
motdir=unique(d);

%strength
coh24=coh==0.24;
coh12=coh==0.12;
coh06=coh==0.06;
kl(coh==0.24)=simP(1);
kl(coh==0.12)=simP(2);
kl(coh==0.06)=simP(3);
m=1:1:360;

%1.Bayesian inference with sensory evidence and cardinal priors
%--------------------------------------------------------------
%get a matrix (360 possible MAPs, 360 possible motion directions) of 
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4 
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.
%get likelihood of data in condition with coherence 24
[~,likeMAP11]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior1,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP12]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior2,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP13]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior3,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP14]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior4,0,0,kcardinal,priorShape,TheModel,varargin);

%get likelihood of data in condition with coherence 12
[~,likeMAP21]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior1,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP22]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior2,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP23]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior3,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP24]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior4,0,0,kcardinal,priorShape,TheModel,varargin);

%get likelihood of data in condition with coherence 6
[~,likeMAP31]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior1,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP32]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior2,0,0,kcardinal,priorShape,TheModel,varargin);
[~,likeMAP33]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior3,0,0,kcardinal,priorShape,TheModel,varargin);
[MAP,likeMAP34]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior4,0,0,kcardinal,priorShape,TheModel,varargin);

%get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials) 
%for possible values of upo (rows) for each trial (column)
%Probability of getting a given estimate given Bayesian inference
%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials) 
%for possible values of upo (rows) for each trial (column)
PupoGivenBIwithEvandCard=nan(numel(MAP),numel(d));
for i=1:numel(motdir)
    
    %get displayed motion direction for this trial
    thisd=d==motdir(i);
    
    %get likelihood of data in condition with learnt prior 80
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior1)=likeMAP11(:,motdir(i(ones(sum(thisd&coh24&Prior1),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior1)=likeMAP21(:,motdir(i(ones(sum(thisd&coh12&Prior1),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior1)=likeMAP31(:,motdir(i(ones(sum(thisd&coh06&Prior1),1))));
    
    %get likelihood of data in condition with learnt prior 40
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior2)=likeMAP12(:,motdir(i(ones(sum(thisd&coh24&Prior2),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior2)=likeMAP22(:,motdir(i(ones(sum(thisd&coh12&Prior2),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior2)=likeMAP32(:,motdir(i(ones(sum(thisd&coh06&Prior2),1))));

    %get likelihood of data in condition with learnt prior 20
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior3)=likeMAP13(:,motdir(i(ones(sum(thisd&coh24&Prior3),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior3)=likeMAP23(:,motdir(i(ones(sum(thisd&coh12&Prior3),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior3)=likeMAP33(:,motdir(i(ones(sum(thisd&coh06&Prior3),1))));
    
    %get likelihood of data in condition with learnt prior 10
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior4)=likeMAP14(:,motdir(i(ones(sum(thisd&coh24&Prior4),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior4)=likeMAP24(:,motdir(i(ones(sum(thisd&coh12&Prior4),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior4)=likeMAP34(:,motdir(i(ones(sum(thisd&coh06&Prior4),1))));
end

%scale to probability.
Z_=sum(PupoGivenBIwithEvandCard);
Z=Z_(ones(size(PupoGivenBIwithEvandCard,1),1),:);
PupoGivenBIwithEvandCard=PupoGivenBIwithEvandCard./Z;


%1.Bayesian inference with learnt prior and cardinal priors
%----------------------------------------------------------
[~,likeMAP11]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP12]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP13]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP14]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin);

%get likelihood of data in condition with coherence 12
[~,likeMAP21]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP22]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP23]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP24]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin);

%get likelihood of data in condition with coherence 6
[~,likeMAP31]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP32]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin);
[~,likeMAP33]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin);
[MAP,likeMAP34]=SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin);

%Probability of getting a given estimate
%now get matrix 'PupoGivenBIwithlearntAndCard' of likelihood values 
%(upos=1:1:360,trials) for possible values of upo (rows) for each trial 
%(column)
PupoGivenBIwithlearntAndCard=nan(numel(MAP),numel(d));
for i=1:numel(motdir)
    
    %get displayed motion direction for this trial
    thisd=d==motdir(i);
    
    %get likelihood of data in condition with learnt prior 80
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior1)=likeMAP11(:,motdir(i(ones(sum(thisd&coh24&Prior1),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior1)=likeMAP21(:,motdir(i(ones(sum(thisd&coh12&Prior1),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior1)=likeMAP31(:,motdir(i(ones(sum(thisd&coh06&Prior1),1))));
    
    %get likelihood of data in condition with learnt prior 40
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior2)=likeMAP12(:,motdir(i(ones(sum(thisd&coh24&Prior2),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior2)=likeMAP22(:,motdir(i(ones(sum(thisd&coh12&Prior2),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior2)=likeMAP32(:,motdir(i(ones(sum(thisd&coh06&Prior2),1))));

    %get likelihood of data in condition with learnt prior 20
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior3)=likeMAP13(:,motdir(i(ones(sum(thisd&coh24&Prior3),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior3)=likeMAP23(:,motdir(i(ones(sum(thisd&coh12&Prior3),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior3)=likeMAP33(:,motdir(i(ones(sum(thisd&coh06&Prior3),1))));
    
    %get likelihood of data in condition with learnt prior 10
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior4)=likeMAP14(:,motdir(i(ones(sum(thisd&coh24&Prior4),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior4)=likeMAP24(:,motdir(i(ones(sum(thisd&coh12&Prior4),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior4)=likeMAP34(:,motdir(i(ones(sum(thisd&coh06&Prior4),1))));
end

%scale to probability.
Z_=sum(PupoGivenBIwithlearntAndCard);
Z=Z_(ones(size(PupoGivenBIwithlearntAndCard,1),1),:);
PupoGivenBIwithlearntAndCard=PupoGivenBIwithlearntAndCard./Z;

%get probabilities of the range of percepts "upo" given random estimation
%------------------------------------------------------------------------
PupoGivenRand=ones(360,numel(d))/360;

%calculate the fraction of trials that is controlled by each process 
%We use a divisive normalization competition rule between the likelihood, 
%the learnt prior and the cardinal prior. Basically the stronger the 
%representation and the
%most often it will be chosen. We use the strength of the von Mises
%composing the cardinal prior and not some measure of the overall strength
%of the cardinal prior. But both are correlated. For fixed modes, if the
%strength of the von Mises (one value here) increases, so does the strength
%of the overall distribution. I think...
%The probability of random estimation is fixed.
weightPriorCardmn = klearnt./(klearnt+kl);
weightLlhCardmn = kl./(klearnt+kl);

%scale the mixing weights to probabilities (all sum to 1)
%this is not equal to 1. We want to make it equal to 1.
sumP = weightPriorCardmn + weightLlhCardmn + Prandom;
PpriorCardmn = 1 - (Prandom + weightLlhCardmn)./sumP;
PllhCardmn = 1 - (Prandom + weightPriorCardmn)./sumP;
Prandnew = unique(1 - (weightLlhCardmn + weightPriorCardmn)./sumP);

%take just one value of Prandnew because sometimes because of numerical
%instability there is slight variation in the value of Prandnew which
%produces more than one value. Mathematically there should be only one
%value because priors and llh strength have been scaled some that they sum
%to 1 and Prandom is a constant. Thus Prandnew must be one value.
Prandnew=Prandnew(1);

%repeat P(choose prior), matrix values, of each trial (columns) for each 
%possible data value (rows).
numm=numel(m);
PpriorCardmnall=PpriorCardmn(ones(numm,1),:);
PllhCardmnall=PllhCardmn(ones(numm,1),:);

%Either one chooses evidence (llh mean) or he chooses the prior mean. 
%P(choose evidence|motion) = P(obs evidence & choose evidence/motion and prior)
%P(choose prior mean|motion) = P(choose prior mean/motion and prior)
%Thus equal llh and prior strengths do not produces a bimodal estimate
%density with two peaks with same amplitude. It would be the case if either
%subject chose the motion direction (mean of sensory density) or the prior
%mean(Competition between prior and sensory density). Here it's competition
%between likelihood mean and prior mean. An alternative model with make
%motion direction (mean of evidence density) and prior mean compete instead
%of the mean of likelihood (sampled at each trial from evidence density)
%and the prior mean.
PupoGivenModel=PupoGivenBIwithEvandCard.*PllhCardmnall + PupoGivenBIwithlearntAndCard.*PpriorCardmnall + PupoGivenRand.*Prandnew;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel))==1
    keyboard
    fprintf('%s \n','Something s wrong. PupoGivenModel are probabilties and should sum to 1')
end

%convolve with motor noise
%--------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo=0:1:359;
Pmot=vmPdfs(upo,0,km,'norm');
Pmotforconv=Pmot(:,ones(1,numel(d)));
PestimateGivenModel=SLcircConv(PupoGivenModel,Pmotforconv);

%set upo to initial values any case we use it later
upo=1:1:360;

%trial-predictions (sample estimate density). Variability in estimate
%should reflect variability in measurement density and motor noise
for i=1:length(d)
    pred(i)=randsample(1:1:360,1,'true',PestimateGivenModel(:,i));    
end

%draw sampled predictions
%------------------------
%draw predictions and data
%case von Mises Prior
%--------------------
if sum(strcmp(varargin,'vonMisesPrior')) 
    drawCircStat(pred,d,coh,pstd,[]);
    drawCircStat(data,d,coh,pstd,[]);
end

%case bimodal Prior
%--------------------
if sum(strcmp(varargin,'bimodalPrior'))
    drawCircStat(pred,d,coh,priorModes(:,2)-priorModes(:,1),[]);
    drawCircStat(data,d,coh,priorModes(:,2)-priorModes(:,1),[]);
end

%get loglikelihood of data
%-------------------------
%single trial's measurement, its position(row) for each trial(col) and its
%probability (also maxlikelihood of trial's data). Checked many times. It
%works.
%make sure sub2ind inputs are the same size
if sum(size(data)~=size([1:1:numel(d)]'))==2
    data=data';
end
idx=sub2ind(size(PestimateGivenModel),data,[1:1:numel(d)]');
PdataGivenModel=PestimateGivenModel(idx);

%We use log likelihood because likelihood is so small that matlab cannot
%encode it properly (numerical unstability). We can use single trials log
%likelihood to calculate AIC in the conditions that maximize differences in
%predictions of two models.
logL_pertrial=log(PdataGivenModel);

%We use -sum(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
negSumlogL=-sum(log(PdataGivenModel));

%We use E(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
negElogL=-nanmean(log(PdataGivenModel));

%Look at fitting. It is 3X faster without drawing.
ti=toc;
fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
    negElogL,kl1,kl2,kl3,klearnt1,klearnt2,klearnt3,klearnt4,kcardinal,Prandom,km,ti)


%draw model's internal representations and predictions 
%-----------------------------------------------------
%for now
if sum(strcmp(varargin,'vonMisesPrior'))==1
    SLdrawModelRepresentations(coh,pstd,priorModes,d,simP,PestimateGivenModel,TheModel,varargin)
end