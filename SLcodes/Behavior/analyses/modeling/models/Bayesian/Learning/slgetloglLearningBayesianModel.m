
% slgetloglLearningBayesianModel.m
%
%
% author: steeve laquitaine
%   date: 141122 updated 150203
%purpose: get the log likelihood that Bayesian Models produce observed data
%
%reference: Stocker and Simoncelli, 2006, NN
%
%
%usage:
%
%       [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP] = slgetloglLearningBayesianModel(data,displ,coh,pstd,fitP,priorModes,TheModel,'bimodalPrior','MAPReadout')
%
%
%varargin options
%
%                   'FitDisplayOn': plot data and model predictions
%
%Predictions:
%
%Difference in models predictions is maximized when:
%- motion is weak enough to exploit prior
%- prior is weak enough to produce direction distant from the modes but not
%strong enough to compete with motion.



%logl
function [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP,AIC,PdataGivenModel,duration] = slgetloglLearningBayesianModel(data,...
    displ,...
    coh,...
    pstd,...
    fitP,...
    priorModes,...
    trials,...
    TheModel,varargin)

tic

%Data from each trial can be explained by three Bayesian model hypotheses:
%
%1) Bayesian inference (MAP readout)
%----------------------------------
%     1) measurement mi (=llh mean: ul) is drawn from a von Mises measurement
%     density (v(mi,km)) in 1-"Pp"+"Prandom" fraction of trials.
%
%     2) At each trial likelihood peaks at the sampled measurement and is
%     combined with a prior (assumed von Mises for unimodal priors or
%     mixture of von Mises for bimodal priors with the strength a fit
%     parameter).
%     This produce a distribution of estimates for each motion direction.
%
%     3) Random estimation
%     Sometimes, in 'Prandom' fractions of trials, subjects randomly sample
%     a direction estimate from a uniform distribution when he doesn't know
%     what to choose. So the estimate distribution becomes a sum
%     of the Bayesian estimate distribution and the uniform distribution
%     respectively weighted by 1-Prandom and Prandom.
%
%     4) motor noise ~v(0,kmo) is added to the estimate (let's call it percept)
%     of those two processes by convolving motor noise with the previously
%     generated estimate distribution.
%
%     5) prior learning rate parameter tau
%
%Thus, the maximum likelihood of observing each trial data assuming this
%model is true is:
%
%..........................................................................
%  p(estimate|Bayes)*(1-"Pp"+"Prandom") + p(estimate|random)*Prandom
%..........................................................................
%
%convolved with motor noise (note: its is because convolution
%is distributive).
%
%time
ticLogL = tic;
ulall = 1:1:360;

%check number of parameters
checkNumberOfParams(fitP,12)

%getParameters
[kl1,kl2,kl3,kcardinal,Prandom,km,tau] = getParameters(fitP);

%set parameter boundaries
%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time. Constrain other fraction 
%parameters between 0 and 1, and other parameters as>0.
[negLogl,Logl_pertrial,skip] = setParamBoundaries(Prandom,fitP);
if skip==1; return; end

%set LLHs
[numLLH, TailLLH] = setLLHs(coh);

%set priors
[SetOfpriorModes,numPriors,TailPrior,kpriors,PstdSets,Prior,modesPrior] = setPriors(priorModes,pstd,numLLH,fitP);

%1.Bayesian inference
%get a matrix (360 possible MAPs, 360 possible motion directions) of
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.

%We only calculate
%likelihood of the data for the actually displayed motion direction and
%not the full range of 360 motion directions.

%To adjust the code to fit combined priors, we have to change way we
%fit bimodal prior experiment. Instead fitting 4 disting parameters to
%estimate the srengths of the 20 experimental von mises composing the
%4 bimodal priors. We now use only 1 fit parameter for the 4 (corresponding
%to kpriors3 in the von Mises prior case).

%In the case of bimodal prior, the fitting procedure uses only Klearnt3
%the 6th fit parameter to estimate the four prior strength (whatever what
%we input in the 3 other prior strength parameters).

%In the case of von Mises prior, the four kpriors1 to kpriors 3 parameters
%are used to fit the data.

%In the case of combined prior, the three kpriors1, kpriors2, kpriors4
%parameters are used to fit vonMises data, and kpriors3 is used to fit
%both for vonMises
numTrial = nan(1,numPriors);
for i = 1 : numPriors
    
    %get # of trial in prior block
    numTrial(i) = max(trials(pstd==PstdSets(i)));
    thisPrior = pstd==PstdSets(i);
    thisKp = kpriors(i);
        
    %parfor t = 1 : numTrial(i)
    parfor t = 1 : numTrial(i)
        
        %calculate learning factor (trial scaling of kp)
        kpriors_t = (1-exp(-(t-1)/tau))*thisKp;
        
        %find unique motion directions for trial t in prior block i
        motdir = unique(displ(thisPrior & trials==t));
        
        %coh 24,12,6
        [MAP(:,t),likeMAP1(:,:,i,t)] = slLearningBayesLookupTable(ulall,motdir,kl1,TailLLH,modesPrior(i,:),kpriors_t,kcardinal,TailPrior,TheModel);
        [~  ,likeMAP2(:,:,i,t)] = slLearningBayesLookupTable(ulall,motdir,kl2,TailLLH,modesPrior(i,:),kpriors_t,kcardinal,TailPrior,TheModel);
        [~  ,likeMAP3(:,:,i,t)] = slLearningBayesLookupTable(ulall,motdir,kl3,TailLLH,modesPrior(i,:),kpriors_t,kcardinal,TailPrior,TheModel);
    end
end

%LLH
coh24 = coh==0.24;
coh12 = coh==0.12;
coh06 = coh==0.06;

%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials)
%for possible values of upo (rows) for each trial (column)
PupoGivenBI = nan(numel(MAP(:,1)),numel(data));

for j = 1 : numPriors
    for t = 1 : numTrial(j)
        
        thisCol24 = Prior(:,j) & trials==t & coh24;
        motdir24 = displ(thisCol24);        
        thisCol12 = Prior(:,j) & trials==t & coh12;
        motdir12 = displ(thisCol12);        
        thisCol06 = Prior(:,j) & trials==t & coh06;
        motdir06 = displ(thisCol06);
        
        %coh 24,12,6
        PupoGivenBI(:,thisCol24) = likeMAP1(:,motdir24,j,t);
        PupoGivenBI(:,thisCol12) = likeMAP2(:,motdir12,j,t);
        PupoGivenBI(:,thisCol06) = likeMAP3(:,motdir06,j,t);
    end
end

%scale to probability.
%This is to make sure but it should already scale to probability and sum to
%1.
Z_ = sum(PupoGivenBI);
Z =  Z_(ones(size(PupoGivenBI,1),1),:);
PupoGivenBI =  PupoGivenBI./Z;

%probabilities of percepts "upo" given random estimation
%-------------------------------------------------------
PupoGivenRand =  ones(360,numel(displ))/360;

%calculate probability of percepts "upo" given the model
PBI =  1 - Prandom;
PupoGivenModel =  PupoGivenBI*PBI + PupoGivenRand*Prandom;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel))==1
    keyboard
    fprintf('%s \n','Something wrong. PupoGivenModel are probabilties and should sum to 1')
end

%convolve with motor noise
%-------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo = 0:1:359;
Pmot = vmPdfs(upo,0,km,'norm');
Pmotforconv = Pmot(:,ones(1,numel(displ)));
PestimateGivenModel = SLcircConv(PupoGivenModel,Pmotforconv);

%check that probability of estimates Given Model are positive values.
%circular convolution sometimes produces negative values very close to zero 
%(order of -10^-18). Those values produce infinite -log likelihood which are
%only need a single error in estimation to be rejected during model fitting 
%(one trial has +inf -log likelihood to be predicted by the model).
%This is obviously too conservative. We need a minimal non-zero lapse rate. 
%Prandom that allows for errors in estimation. So we add 10^320 the minimal
%natural number available in matlab. This means that every time an estimate
%that is never produced by the model without lapse rate is encountered
%-loglikelihood increases by -log(10^-320) = 737 and is thus less likely to
%be a good model (lowest -logLLH). But the model is not rejected altogether
%(it would be -log(0) = inf). In the end models that cannot account for
%error in estimates are more likely to be rejected than models who can.
if ~isempty(find(PestimateGivenModel<=0))    
    PestimateGivenModel(PestimateGivenModel<=0) = 10^-320;
    fprintf('(SLgetLoglBayesianModel) P(estimates|model) were sometimes negative but very close to zero < 10^-10 and were thus floored at the lowest positive numerical value 10^-320. \n')
end

%set upo to initial values any case we use it later
upo = 1 : 1 : 360;

%normalize to probability; If we don't normalize, more random
%choice always increase probability of observing data causing larger
%probability of random choice to prevail. To avoid that we need to normalize
%to probabilities that sum to 1. It also makes intuitive sense to deal
%with probabilities.
%Here too PestimateGivenModel integrated over estimates should already sum to
%1. But we still scale to be absolutely sure.
Z_ = sum(PestimateGivenModel);
Z = Z_(ones(size(PestimateGivenModel,1),1),:);
PestimateGivenModel = PestimateGivenModel./Z;


%case we want the log likelihood of the data
%-------------------------------------------
%other case are when we just want estimates distributions prediction given
%model parameters.
data(data==0) = 360;
if isempty(data)==0
    
    %single trial's measurement, its position(row) for each trial(col) and its
    %probability (also maxlikelihood of trial's data). Checked many times. It
    %works.
    %make sure sub2ind inputs are the same size
    if sum(size(data)~=size([1:1:numel(displ)]'))==2
        data = data';
    end
    idx = sub2ind(size(PestimateGivenModel),data,[1:1:numel(displ)]');
    PdataGivenModel = PestimateGivenModel(idx);
    
    %check that likelihood are positive 
    if ~isempty(find(PdataGivenModel<=0))        
        fprintf('(SLgetLoglBayesianModel) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        error('(SLgetLoglBayesianModel) !!! WARNING !!! likelihood of data is negative. Something is wrong !')
        fprintf('(SLgetLoglBayesianModel) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')        
    end
    
    %We use log likelihood because likelihood is so small that matlab cannot
    %encode it properly (numerical unstability). We can use single trials log
    %likelihood to calculate AIC in the conditions that maximize differences in
    %predictions of two models.
    Logl_pertrial = log(PdataGivenModel);
    
    %check that logl are real numbers
    if ~isreal(PdataGivenModel)        
        fprintf('(SLgetLoglBayesianModel) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        error('(SLgetLoglBayesianModel) !!! WARNING !!! log likelihood of data are not real numbers. Something is wrong !')
        fprintf('(SLgetLoglBayesianModel) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')        
    end
    
    %We use -sum(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.                
    negLogl = round(-sum(Logl_pertrial));        
    
    %AIC
    nfitP = sum(~isnan(fitP));
    AIC = 2 * (nfitP - sum(Logl_pertrial));
    
    %We use -sum(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    %negSumlogL = -sum(log(PdataGivenModel));
    
    %We use E(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    %negElogL = -nanmean(log(PdataGivenModel));    

    %print status
    duration = toc;
    fp = repmat('%.2f ',1,12);    
    fprintf(['%.2f %s' fp '%s %.2f%s \n'],negLogl,'[',fitP,']',duration,'s')
    
    %display fit
    %Look at fitting. It is 3X faster without drawing.   
%     if sum(strcmp(varargin{:},'FitDisplayOn'))
%         
%         %mean and std
%         %data
%         [meanData,stdData,dataCond] = SLmakeDataMeanAndStd(data,displ,...
%             coh,pstd,priorModes,priorShape);
%         
%         %pred
%         [meanPred, stdPred] = makePred(displ,coh,pstd,priorShape,priorModes,'Mean',[],PestimateGivenModel,MAP);
%         
%         %draw
%         SLdrawModelsPredictionCentered(meanData,stdData,meanPred,stdPred,dataCond,...
%             priorModes,priorShape,'yCentered')
%     end
    
    %display parameters
    ti = toc(ticLogL);    
end


%predictions
function [meanPred, stdPred,skip] = makePred(displ,coh,pstd,priorShape,priorModes,TrialOrMean,output,PestimateGivenModel,MAP)

%case von Mises prior
%--------------------
%conditions
if strcmp(priorShape,'vonMisesPrior')
    [cond,idxCondUniqtrial,~] = SLuniqpair([pstd coh displ]);
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
if strcmp(TrialOrMean,'Trial')    
    %von Mises Prior    
    if strcmp(priorShape,'vonMisesPrior')
        output.TrialPred = nan (size(pstd));
        for i = 1 : numCond
            trialsThiC = pstd==cond(i,1) & coh==cond(i,2) & displ==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
end

function checkNumberOfParams(fitP,nb)

%warning
if length(fitP) < nb
    fprintf('%s \n','(getLogl) Some fit parameters are missing')
    keyboard
end

%get parameters
function [kl1,kl2,kl3,kcardinal,Prandom,km,tau] = getParameters(fitP)

%They should be in this order:
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal',
%'Prand','km','tail','tau'
kl1       = fitP(1);
kl2       = fitP(2);
kl3       = fitP(3);
kcardinal = fitP(8);
Prandom   = fitP(9);
km        = fitP(10);
tau       = fitP(12);

%no cardinal prior
if isnan(kcardinal)
    kcardinal = 0;
end

%set parameters boundaries
function [negLogl,Logl_pertrial,skip] = setParamBoundaries(Prandom,fitP)

negLogl = []; Logl_pertrial = [];

%case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom > 1
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    fprintf('%s \n','(getLogl) Prandom is > 1 - skip')
    return
end
skip=0;
if any(fitP < 0)
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    fprintf('%s \n','(getLogl) a fit parameter was > 1 - skip')
    skip=1;
    return
end
if any(isnan(fitP([1,2,3,4,5,6,7,9,10,12])))
    negLogl = 1e09;
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal or TailPrior is NaN'])
    keyboard
end

%set LLHS
function [numLLH, TailLLH] = setLLHs(coh)

numLLH = numel(unique(coh));
TailLLH = zeros(numLLH,1);

%set priors
function [SetOfpriorModes,numPriors,TailPrior,kpriors,PstdSets,Prior,modesPrior] = setPriors(priorModes,pstd,numLLH,fitP)

SetOfpriorModes = SLuniqpair(priorModes);
numPriors = numel(unique(pstd));
TailPrior = NaN;
for i = 1 : numPriors; kpriors(i) = fitP(numLLH + i); end
PstdSets = sort(unique(pstd),'descend');
Prior = nan(numel(pstd),numPriors);
for i = 1 : numPriors
    Prior(:,i) = pstd==PstdSets(i);
    modesPrior(i) = SetOfpriorModes;
end
modesPrior = SLmakeColumn(modesPrior);

