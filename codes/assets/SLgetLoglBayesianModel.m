
%SLgetLoglBayesianModel.m
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
%       [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP] = SLgetLoglBayesianModel(data,displ,coh,pstd,fitP,'bimodalPrior',priorModes,TheModel,'bimodalPrior','MAPReadout')
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


function [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP,AIC,PdataGivenModel] = SLgetLoglBayesianModel(data,...
    displ,...
    coh,...
    pstd,...
    fitP,...
    priorShape,...
    priorModes,...
    TheModel,...
    varargin)

%time
ticLogL = tic;

%case it is the motion energy model set motion
%strength parameters (sensory likelihood was only computed
%for )
%LLHs strengths are always the numLLHs first
%parameters
vrg = varargin;
numVarg = length(vrg{:});
if slIsInput(vrg{:},'motionenergymodel')
    vrg{:}(numVarg+1) = {'motionStrength'};
    numLLHs = 0;
    kl = [nan nan nan];
else
    vrg{:}(numVarg+1) = {'stimStrength'};
    numLLHs = numel(unique(coh));
    kl = fitP(1:numLLHs);
end
stimStrength = sort(unique(coh),'descend');

%and priors the next parameters
%and 4 always fixed other parameters
%kcardinal, prandom, km and prior tail mixture coefficient
numPriors = numel(unique(pstd));

%(case von Mises prior)
%----------------------
SetOfpriorModes = SLuniqpair(priorModes);
if sum(strcmp(priorShape,'vonMisesPrior'))
    
    %case stimulus and prior strengths conditions
    %are missing in the subset of the dataset and replaced
    %by NaN.
    %remove parameters from the stimulus strength and
    %the prior strength that are NaN from the fit parameters.
    %This assumes the last 3 parameters are always the cardinal
    %prior, prandom and prior tail
    %case a model that is not the motion energy model
    if ~slIsInput(vrg{:},'motionenergymodel')
        fitPLLH_Prior = fitP(1:length(fitP)-4);
        if any(isnan(fitPLLH_Prior))
            fitPLLH_Prior(isnan(fitPLLH_Prior)) = [];
            fitP = [fitPLLH_Prior fitP(length(fitP)-3:end)];
        end
    elseif slIsInput(vrg{:},'motionenergymodel')
        fitP_Prior = fitP(1:numPriors);
        if any(isnan(fitP_Prior))
            fitP_Prior(isnan(fitP_Prior)) = [];
            fitP = [fitP_Prior fitP(end-3:end)];
        end
    end
    
    %scale the parameters to have the same units
    %improve fitting (picking up the scaling parameters
    %assigned to LLHs and priors)
    if slIsInput(vrg{1},'pscaling')
        pscales = vrg{1}{find(strcmp(vrg{:},'pscaling'))+1};
        pscales = [pscales(1:numLLHs) pscales(numLLHs+1:numLLHs+numPriors) pscales(end-3:end)];
        fitP = fitP./pscales;
    end
    
    %Trial matrix of stim strength conditions
    %each column of LLHs corresponds to a row of
    %stimStrengthSets
    LLHs = nan(numel(coh),numLLHs);
    stimStrengthSets = sort(unique(coh),'descend');
    for i = 1 : length(stimStrength)
        LLHs(:,i) = coh==stimStrengthSets(i);
    end
    
    %case fat tail
    if sum(strcmp(vrg{:},'FatTailPrior'))==1
        
        [kpriors,TailPrior] = setupFatTailPrior(fitP,numLLHs,numPriors,vrg);
        
        %(case fat tail priors and llhs)
    elseif sum(strcmp(vrg{:},'FatTailPriorAndLLH'))==1 && sum(strcmp(vrg{:},'FitEachKvmAndTails'))==0
        
        %case tail fixed across priors, Kvm free
        %tails and fitP 4 to ... are priors
        if sum(strcmp(vrg{:},'ChangeWTail'))==0
            TailPrior = fitP(11);
            for i = 1 : numPriors
                kpriors(i) = fitP(3+i);
            end
        end
        
        %case tail free across priors, Kvm fixed
        if sum(strcmp(vrg{:},'ChangeWTail'))==1
            %tails
            %fitP 4 to ... are priors
            %k priors
            for i = 1 : numPriors
                TailPrior(i) = fitP(3+i);
            end
            kpriors = fitP(11);
        end
        
        %case fat tail prior and llh, all tails and Kvm free
    elseif sum(strcmp(vrg{:},'FatTailPriorAndLLH'))==1 && sum(strcmp(vrg{:},'FitEachKvmAndTails'))==1
        
        %tails
        %WARNING: Make sure the fit parameters are in the following
        %order: kl24, kl12, kl6, kpriors1, kpriors2, kpriors3, kpriors4, kcard, Prand, Kmotor,
        %Tail_llh1, Tail_llh2, Tail_llh3, Tail_prior1, Tail_prior2, Tail_prior3, Tail_prior4
        
        %fitP numLLH to ... are priors
        for i = 1 : numPriors
            kpriors(i) = fitP(numLLHs+i);
        end
        
        %fitP 11 to ... are llh tails
        for i = 1 : numLLHs
            TailLLH(i) = fitP(numLLHs + numPriors + 3 + i);
            %case Tail is fixed set TaiLLH = 0 if NaN
            if isnan(TailLLH(i))
                TailLLH(i) = 0;
            end
        end
        
        %fitP 14 to ... are prior tails
        for i = 1 : numPriors
            TailPrior(i) = fitP(numLLHs + numPriors + 3 + numLLHs + i);
            %case Tail is fixed set TailPrior = 0 if NaN
            if isnan(TailPrior(i))
                TailPrior(i) = 0;
            end
        end
        
        %(case von Mises without fat tail)
    elseif sum(strcmp(vrg{:},'FatTailPrior'))==0 && sum(strcmp(vrg{:},'FatTailPriorAndLLH'))==0
        %tails (should be NaN)
        %fitP 4 to ... are priors
        TailPrior = NaN;
        for i = 1 : numPriors
            kpriors(i) = fitP(numLLHs + i);
        end
        TailLLH = zeros(length(stimStrength),1);
    end
    
    %trials
    %fitP 4 to n... are the priors
    PstdSets = sort(unique(pstd),'descend');
    Prior = nan(numel(pstd),numPriors);
    for i = 1 : numPriors
        Prior(:,i) = pstd==PstdSets(i);
        modesPrior(i) = SetOfpriorModes;
    end
    modesPrior = SLmakeColumn(modesPrior);
    
    %fitP 4 to 7 are priors
    %get fitP assigned to the priors present in the dataset.
    kpriors = fitP(numLLHs+1:numLLHs+numPriors);
end

%warning
if length(fitP) < numLLHs+numPriors+4
    fprintf('%s \n','(getLogl) Some fit parameters are missing')
    keyboard
end
ulall = 1 : 1 : 360;

%other fitp always come after the likelihood
%and prior parameters in the same order
kcardinal = fitP(numLLHs+numPriors+1);
Prandom   = fitP(numLLHs+numPriors+2);
km        = fitP(numLLHs+numPriors+3);

%output
negLogl       = [];
Logl_pertrial = [];

%------
%PRIORS
%------
%(case cardinal priors)
%parameter is fixed
if isnan(kcardinal)
    kcardinal = 0;
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
if ~slIsInput(varargin{:},'motionenergymodel') && any(isnan([kl kpriors Prandom km]))
    negLogl = 1e09;
    dbstack
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal or TailPrior is NaN'])
    keyboard
end

%case no fat tail
if sum(isnan(TailPrior))==1
    TailPrior = 0;
end

%warning
if any(TailPrior > 1)
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    fprintf('%s \n','(SLgetLoglBayesianModel) TailPrior is > 1 - skip')
    return
end

%1.Bayesian inference
%get a matrix (360 possible MAPs, 360 possible motion directions) of
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.

%motion directions (not fixed)
%resist to change in the dataset (e.g., if one
%only take a subset of the data)
motdir = unique(displ);

%------------------
%(case MAP readout)
%------------------
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
%both for vonMises and bimodal prior data.
if sum(strcmp(varargin{:},'MAPReadout'))
    for i = 1 : numPriors
        for j = 1 : length(stimStrength)
            %record stimulus strength
            vrg{:}{numVarg+2} = stimStrength(j);
            %compute percept distribution
            [MAP,likeMAP{i}{j}] = SLGirshickBayesLookupTable(ulall,motdir,kl(j),modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,vrg{:});
        end
    end
end

%(case Sampling readout)
%-----------------------
%likelihood of data
if sum(strcmp(varargin{:},'SamplingReadout'))==1
    %case priors have same tail weight
    if numel(TailPrior)==1
        for i = 1 : numPriors
            for j = 1 : length(stimStrength)
                [MAP,likeMAP{i}{j}] = SLBayesSamplingLookupTable(ulall,motdir,kl(j),TailLLH(1),modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,vrg{:});
            end
        end
    end
    if numel(TailPrior)>1
        for i = 1 : numPriors
            for j = 1 : numLLHs
                [MAP,likeMAP{i}{j}] = SLBayesSamplingLookupTable(ulall,motdir,kl(j),TailLLH(1),modesPrior(i,:),kpriors(i),kcardinal,TailPrior(i),priorShape,TheModel,vrg{:});
            end
        end
    end
end

%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials)
%for possible values of upo (rows) for each trial (column)
PupoGivenBI = nan(numel(MAP),numel(data));
for i = 1 : numel(motdir)
    %displayed direction
    thisd = displ==motdir(i);
    for j = 1 : numPriors
        for k = 1 : length(stimStrength)
            PupoGivenBI(:,thisd&LLHs(:,k)&Prior(:,j)) = likeMAP{j}{k}(:,motdir(i(ones(sum(thisd&LLHs(:,k)&Prior(:,j)),1))));
        end
    end
end
clear likeMAP;

%scale to probability.
%This is to make sure but it should already scale to probability and sum to
%1.
Z_ = sum(PupoGivenBI);
Z =  Z_(ones(size(PupoGivenBI,1),1),:);
PupoGivenBI =  PupoGivenBI./Z;

%probabilities of percepts "upo" given random estimation
PupoGivenRand =  ones(360,numel(displ))/360;

%calculate probability of percepts "upo" given the model
PBI =  1 - Prandom;
PupoGivenModel =  PupoGivenBI*PBI + PupoGivenRand*Prandom;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel)) ==1
    keyboard
    fprintf('%s \n','(SLgetLoglBayesianModel) PupoGivenModel should sum to 1')
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
    negLogl = -sum(Logl_pertrial);
    
    %AIC
    nfitP = sum(~isnan(fitP));
    AIC = 2 * (nfitP - sum(Logl_pertrial));
    
    %We use -sum(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    %negSumlogL = -sum(log(PdataGivenModel));
    
    %We use E(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    %negElogL = -nanmean(log(PdataGivenModel));
    
    %Case not both priors combined.
    %Look at fitting. It is 3X faster without drawing.
    
    %display fit
    if sum(strcmp(varargin{:},'FitDisplayOn'))
        
        %mean and std
        %data
        [meanData,stdData,dataCond] = SLmakeDataMeanAndStd(data,displ,...
            coh,pstd,priorModes,priorShape);
        
        %pred
        [meanPred, stdPred] = makePred(displ,coh,pstd,priorShape,priorModes,'Mean',[],PestimateGivenModel,MAP);
        
        %draw
        SLdrawModelsPredictionCentered(meanData,stdData,meanPred,stdPred,dataCond,...
            priorModes,priorShape,'yCentered')
    end
    
    %display parameters
    ti = toc(ticLogL);
    
    if sum(strcmp(varargin{1},'vonMisesPrior')) && sum(strcmp(varargin{1},'bimodalPrior'))
    else
        %case priors has heavy tail and same von mises strengths with different tail weight
        if sum(strcmp(varargin{:},'ChangeWTail'))==0
            if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==0
                fprintf(['%.05f ',repmat(' %.05f  ',1,length(stimStrength)),repmat(' %.05f  ',1,numPriors),'%.05f  %.05f %.05f %.05f %.1f \n'],...
                    negLogl,kl,kpriors,kcardinal,Prandom,km,TailPrior,ti)
            end
        end
        
        %case priors has heavy tail and same von mises strengths with different tail weight
        if sum(strcmp(varargin{:},'ChangeWTail'))==1
            fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.10f  %.2f  \n',...
                negLogl,kl1,kl2,kl3,TailPrior,kcardinal,Prandom,km,kpriors,ti)
        end
        
        %case priors and llh has free fat tails and von mises strengths
        if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==1
            fprintf('%s %.2f  %.2f %.2f %.2f  %.2f  %.2f  %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.10f  %.10f  %.2f  \n',...
                '(SLgetLoglBayesianModel)',negLogl,kl1,kl2,kl3,kpriors,kcardinal,Prandom,km,TailLLH,TailPrior,ti)
        end
    end
end


%predictions
function [meanPred,stdPred] = makePred(displ,coh,pstd,priorShape,priorModes,TrialOrMean,output,PestimateGivenModel,MAP)

%case von Mises prior
%--------------------
%conditions
if strcmp(priorShape,'vonMisesPrior')
    [cond,idxCondUniqtrial,~] = SLuniqpair([pstd coh displ]);
end

%case bimodal prior
%------------------
%conditions
if strcmp(priorShape,'bimodalPrior')
    
    %check all prior conditions are there
    priorCond = priorModes(:,2) - priorModes(:,1);
    numPriorCond = size(SLuniqpair(priorModes),1);
    
    if numPriorCond == numel(unique(priorCond))
        
        [cond,idxCondUniqtrial,~] = SLuniqpair([priorCond coh displ]);
        
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
if strcmp(TrialOrMean,'Trial')
    
    %case von Mises Prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        output.TrialPred = nan (size(pstd));
        for i = 1 : numCond
            trialsThiC = pstd==cond(i,1) & coh==cond(i,2) & displ==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
    
    %case bimodal prior
    %------------------
    if strcmp(priorShape,'bimodalPrior')
        output.TrialPred = nan(size(priorCond));
        for i = 1 : numCond
            trialsThiC = priorCond==cond(i,1) & coh==cond(i,2) & displ==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
end


%setup fait tail prior model
function [kpriors,TailPrior] = setupFatTailPrior(fitP,numLLHs,numPriors,varargin)

%case tail fixed across priors, Kvm free
%get prior strengths and tail
if sum(strcmp(varargin{:},'ChangeWTail'))==0
    kpriors = fitP(numLLHs+1:numLLHs+numPriors);
    TailPrior = fitP(numLLHs+numPriors+4);
end

%case tail free across priors, Kvm fixed
%tails; fitP 4 to ... are priors; k priors
if sum(strcmp(varargin{:},'ChangeWTail'))==1
    for i = 1 : numPriors
        TailPrior(i) = fitP(3+i);
    end
    kpriors = fitP(11);
end