
%SLmakePredictionsCompetitionModel.m
%
% author: steeve laquitaine
%   date: 141122 updated 150110
%purpose: make Competition model motion estimates predictions.
%
%usage:
%
%       SLmakePredictionsCompetitionModel(d,...
%           coh,pstd,fitP,priorShape,priorModes,[],output,'vonMisesPrior')
%

function [meanPred,stdPred,cond,PestimateGivenModelUniq,MAP,output] = SLmakePredictionsCompetitionModel(d,...
    coh,pstd,fitP,priorShape,priorModes,TrialOrMean,output,varargin)

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

%get fit parameters. They should be in this order.
%(von Mises prior)
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal','Prand','km'
%(bimodal prior)
%'coh24','coh12','coh6','[145 305]','[165 285]','[185 265]','[205 245]',...
%'kcardinal','Prand','km'
tic

%warning
if length(fitP) < 10
    fprintf('%s \n','(getLogl) Some fit parameters are missing')
    keyboard
end

kl1 = fitP(1);
kl2 = fitP(2);
kl3 = fitP(3);
% klearnt1 = fitP(4);
% klearnt2 = fitP(5);
% klearnt3 = fitP(6);
% klearnt4 = fitP(7);
kcardinal = fitP(8);
Prandom = fitP(9);
km = fitP(10);

ulall = 1:1:360;

%initialize
output.TrialPred = [];

%--------------
%SET THE PRIORS
%--------------
%
%(case cardinal prior)
%---------------------
%case we don't want to fit the cardinal prior.
if isnan(kcardinal)
    TheModel = 'withoutCardinal';
    kcardinal=0;
else
    %case we don't want to fit the cardinal prior.
    TheModel='withCardinal';
end


%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time.
%Constrain other fraction parameters between 0 and 1, and other parameters
%as>0.

%case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom > 1
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    fprintf('%s \n','(getLogl) Prandom is > 1 - skip')
    return
end
if any(fitP < 0)
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    fprintf('%s \n','(getLogl) a fit parameter was > 1 - skip')
    return
end
if any(isnan(fitP([1,2,3,4,5,6,7,9,10])))
    negLogl = 1e09;
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal is NaN'])
    keyboard
end

%get prior modes
SetOfpriorModes = SLuniqpair(priorModes);

%(case von Mises prior)
%----------------------
SetOfpriorModes = SLuniqpair(priorModes);
if sum(strcmp(priorShape,'vonMisesPrior')) == 1
    
    %fitP 4 to ... are priors
    numPriors = numel(unique(pstd));
    for i = 1 : numPriors
        klearnt(i) = fitP(3+i);
    end
    
    %trials
    %fitP 4 to n... are the priors 
    PstdSets = sort(unique(pstd),'descend');
    modesPrior = nan(numPriors,1);
    Prior = nan(numel(pstd),numPriors);
    for i = 1 : numPriors
        modesPrior(i) = SetOfpriorModes;
        Prior(:,i) = pstd==PstdSets(i);
    end   
end

%(case bimodal prior)
%--------------------
if sum(strcmp(priorShape,'bimodalPrior')) == 1
    
    %modes and trials
    numPriors = size(SetOfpriorModes,2);
    for i = 1 : numPriors
        modesPrior(i,:) = SetOfpriorModes(i,:);
        Prior(:,i) = priorModes(:,1)==modesPrior(i,1) & priorModes(:,2)==modesPrior(i,2);
    end  
    
    %fitP 4 to ... are priors
    for i = 1 : numPriors
        klearnt(i) = fitP(3+i);
    end 
end


%-----------
%SET THE LLH
%-----------
%get coherences and strength of evidence
%mean
motdir = unique(d);

%strength
coh24=coh==0.24;
coh12=coh==0.12;
coh06=coh==0.06;
m=1:1:360;

%------------------------------------------------------------
%Bayesian inference with sensory evidence and cardinal priors
%------------------------------------------------------------
%get a matrix (360 possible MAPs, 360 possible motion directions) of 
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4 
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.
% We only calculate likelihood of the
%data for the actually displayed motion direction and not the full range of
%360 motion directions. It is faster.
%get likelihood of data in condition with coherence 24
for i = 1 : numPriors
    
    %coh 24,12,6
    [MAP,likeMAP1{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior(i,:),0,kcardinal,0,priorShape,TheModel,varargin{:});
    [~,likeMAP2{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior(i,:),0,kcardinal,0,priorShape,TheModel,varargin{:});
    [~,likeMAP3{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior(i,:),0,kcardinal,0,priorShape,TheModel,varargin{:});
end


%now get matrix of likelihood (upos=1:1:360,trials)
%for possible data (rows) for each trial (column)
PupoGivenBIwithEvandCard = nan(numel(MAP),numel(d));
for i = 1 : numel(motdir)
    
    %displayed motion direction for this trial
    thisd = d==motdir(i);
    
    for j = 1 : numPriors
        
        %coh 24,12,6
        PupoGivenBIwithEvandCard(:,thisd&coh24&Prior(:,j)) = likeMAP1{j}(:,motdir(i(ones(sum(thisd&coh24&Prior(:,j)),1))));        
        PupoGivenBIwithEvandCard(:,thisd&coh12&Prior(:,j)) = likeMAP2{j}(:,motdir(i(ones(sum(thisd&coh12&Prior(:,j)),1))));
        PupoGivenBIwithEvandCard(:,thisd&coh06&Prior(:,j)) = likeMAP3{j}(:,motdir(i(ones(sum(thisd&coh06&Prior(:,j)),1))));
    end
end

%probability.
Z_ = sum(PupoGivenBIwithEvandCard);
Z = Z_(ones(size(PupoGivenBIwithEvandCard,1),1),:);
PupoGivenBIwithEvandCard = PupoGivenBIwithEvandCard./Z;

%--------------------------------------------------------
%Bayesian inference with learnt prior and cardinal priors
%--------------------------------------------------------
for i = 1 : numPriors
    
    %coh 24,12,6
    [MAP,likeMAP{i}] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior(i,:),klearnt(i),kcardinal,0,priorShape,TheModel,varargin{:});
end

%Probability of getting a given estimate
%now get matrix 'PupoGivenBIwithlearntAndCard' of likelihood values 
%(upos=1:1:360,trials) for possible values of upo (rows) for each trial 
%(column)
PupoGivenBIwithlearntAndCard = nan(numel(MAP),numel(d));
for i = 1 : numel(motdir)
    
    %displayed motion direction for this trial
    thisd = d==motdir(i);
    
    for j = 1 : numPriors
        PupoGivenBIwithlearntAndCard(:,thisd&Prior(:,j)) = likeMAP{j}(:,motdir(i(ones(sum(thisd&Prior(:,j)),1))));        
    end
end

%probability.
Z_ = sum(PupoGivenBIwithlearntAndCard);
Z = Z_(ones(size(PupoGivenBIwithlearntAndCard,1),1),:);
PupoGivenBIwithlearntAndCard = PupoGivenBIwithlearntAndCard./Z;

%get probabilities of the range of percepts "upo" given random estimation
%------------------------------------------------------------------------
PupoGivenRand = ones(360,numel(d))/360;


%strengths at each trials
%llh
kl(coh==0.24) = fitP(1);
kl(coh==0.12) = fitP(2);
kl(coh==0.06) = fitP(3);

%priors
PstdSets = sort(unique(pstd),'descend');
klearntTrials = nan(numel(d),1);
for i = 1: numPriors
    klearntTrials(pstd==PstdSets(i)) = klearnt(i);
end
klearntTrials = klearntTrials';

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

%case von Mises Prior
if sum(strcmp(priorShape,'vonMisesPrior')) == 1   
    weightPriorCardmn = klearntTrials./(klearntTrials + kl);
    weightLlhCardmn   = kl./(klearntTrials + kl);
end

%case bimodal prior
%each prior mode plays with its own strength in the competition
if sum(strcmp(priorShape,'bimodalPrior')) == 1    
    %weightPriorCardmn = klearntTrials./(klearntTrials + klearntTrials + kl);
    %weightLlhCardmn   = kl./(klearntTrials + klearntTrials + kl);
    weightPriorCardmn = klearntTrials./(klearntTrials + kl);
    weightLlhCardmn   = kl./(klearntTrials + kl);
end

%scale the mixing weights to probabilities (all sum to 1)
%this is not equal to 1. We want to make it equal to 1.
sumP = weightPriorCardmn + weightLlhCardmn + Prandom;
PpriorCardmn  = 1 - (Prandom + weightLlhCardmn)./sumP;
PllhCardmn    = 1 - (Prandom + weightPriorCardmn)./sumP;
Prandnew      = unique(1 - (weightLlhCardmn + weightPriorCardmn)./sumP);

%take just one value of Prandnew because sometimes because of numerical
%instability there is slight variation in the value of Prandnew which
%produces more than one value. Mathematically there should be only one
%value because priors and llh strength have been scaled some that they sum
%to 1 and Prandom is a constant. Thus Prandnew must be one value.
Prandnew = Prandnew(1);

%repeat P(choose prior), matrix values, of each trial (columns) for each 
%possible data value (rows).
numm = numel(m);
PpriorCardmnall = PpriorCardmn(ones(numm,1),:);
PllhCardmnall = PllhCardmn(ones(numm,1),:);

%We assume here a particular model of competition (llh mean vs prior mean).
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
PupoGivenModel = PupoGivenBIwithEvandCard.*PllhCardmnall + PupoGivenBIwithlearntAndCard.*PpriorCardmnall + PupoGivenRand.*Prandnew;

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
upo = 0:1:359;
Pmot = vmPdfs(upo,0,km,'norm');
Pmotforconv = Pmot(:,ones(1,numel(d)));
PestimateGivenModel = SLcircConv(PupoGivenModel,Pmotforconv);

%get read of extremely small negative numbers that occurs because of 
%instability numerical instability. Set them to zero.
PestimateGivenModel(PestimateGivenModel<0) = 0;

%get predictions densities sorted per condition
%We extract the predicted densities of estimates for each single 
%condition "cond" of the experiment.
%cond(idxCond(thistrial),:) tells the condition (coh pstd d) for each 
%trial SLuniqpair by default gives the index of the first trial with the
%condition in "cond".
%"PestimateGivenModelUniq"'s rows are possible estimates and columns are
%experimental condition in rows of "cond".

%case von Mises prior
%--------------------
%conditions
if strcmp(priorShape,'vonMisesPrior')
    [cond,idxCondUniqtrial,~] = SLuniqpair([pstd coh d]);
end

%case bimodal prior
%------------------
%conditions
if strcmp(priorShape,'bimodalPrior')
    
    %check all prior conditions are there
    priorCond = priorModes(:,2) - priorModes(:,1);
    numPriorCond = size(SLuniqpair(priorModes),1);
    if numPriorCond == numel(unique(priorCond))
        [cond,idxCondUniqtrial,~] = SLuniqpair([priorCond coh d]);
    end
end
PestimateGivenModelUniq = PestimateGivenModel(:,idxCondUniqtrial);

%sort everything to match the data sorting (ascending order)
[cond,PosSorted] = sortrows(cond,[-2 -1]);
PestimateGivenModelUniq = PestimateGivenModelUniq(:,PosSorted);


%predictions about estimate mean and std (circular mean and std) for
%each task condition (columns)
%-------------------------------------------------------------------
%Use circular statistics
numCond = size(PestimateGivenModelUniq,2);
meanPred = nan(numCond,1);
stdPred = nan(numCond ,1);
for i = 1 : numCond 
    data = SLcircWeightedMeanStd(MAP, PestimateGivenModelUniq(:,i));
    meanPred(i)= data.deg.mean;
    stdPred(i) = data.deg.std;
end

%case trial-predictions
%----------------------
%trial-predictions are ordered same as data, pstd, coh and d.
if strcmp(TrialOrMean,'Trial')
    output.TrialPred = nan(size(pstd));
    for i = 1 : numCond
        trialsThiC = pstd==cond(i,1) & coh==cond(i,2) & d==cond(i,3);
        output.TrialPred(trialsThiC) = meanPred(i);
    end
end