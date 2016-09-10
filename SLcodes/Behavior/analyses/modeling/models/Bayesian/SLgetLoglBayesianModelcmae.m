
%SLgetLoglBayesianModelcmae.m
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
%       [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP] = SLgetLoglBayesianModelcmae(data,displ,coh,pstd,fitP,'bimodalPrior',priorModes,TheModel,'bimodalPrior','MAPReadout')
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
function [negLogl,fitP,Logl_pertrial,PestimateGivenModel,MAP,AIC,PdataGivenModel] = SLgetLoglBayesianModelcmae(fitP,data,...
    feature,...
    stimStrength,...
    pstd,...
    priorShape,...
    priorModes,...
    TheModel,...
    varargin)

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
%
%2) Bayesian inference (Sampling readout)
%----------------------------------------
%Model structure is the same but the estimate for each featureayed motion
%direction is sampled from the posterior over repeated trials of each
%motion direction.
%
%
%3) Bayesian inference (Flat tail prior)
%---------------------------------------
%Model structure is the same but the prior is now the mixture of a uniform
%distribution and a von Mises, for unimodal priors, or a mixture of von
%mises, for the bimodal priors.
%We want to test if this model can explain the bimodal estimate
%distribution observed in the data which can not be explained by the
%conventional Bayesian model.

%4) Bayesian inference (Flat tail prior can change)
%--------------------------------------------------

%time
ticLogL = tic;

%rearrange fitP if no cardinal prior
if strcmp(TheModel,'withoutCardinal')
    fitP(end+1)= NaN;
    fitP(8:10) = [NaN;fitP(8:9)]; 
end

%check parameters
if length(fitP) < 10
    fprintf('%s \n','(getLogl) Some fit parameters are missing')
    keyboard
end

%add 11 parameter (tail)
if length(fitP)==10
    fitP(end+1) = NaN;
end

%fix tail weight
if fitP(11)<0
    fitP(11) = NaN;
end
ulall = 1 : 1 : 360;

%fit parameters
%in this order:'c24','c12','c6','p80','p40','p20','p10','kcard',
%'Prand','km'
kl1 = fitP(1);        %coh24
kl2 = fitP(2);        %coh12
kl3 = fitP(3);        %coh06
kcardinal = fitP(8);  %cardinal
Prandom   = fitP(9);  %rand
km        = fitP(10); %motor

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
if any(isnan(fitP([1,2,3,4,5,6,7,9,10])))
    negLogl = 1e09;
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal or TailPrior is NaN'])
    keyboard
end

%(case von Mises prior)
%----------------------
SetOfpriorModes = SLuniqpair(priorModes);
if sum(strcmp(priorShape,'vonMisesPrior'))
    
    %number of priors and LLh
    numPriors = numel(unique(pstd));
    numLLH = numel(unique(stimStrength));

    %case fat tail
    if sum(strcmp(varargin{:},'FatTailPrior'))==1
        
        %case tail fixed across priors, Kvm free
        if sum(strcmp(varargin{:},'ChangeWTail'))==0
            
            %tails
            TailPrior = fitP(11);
            
            %fitP 4 to ... are priors
            for i = 1 : numPriors
                kpriors(i) = fitP(3+i);
            end
        end
        
        %case tail free across priors, Kvm fixed
        if sum(strcmp(varargin{:},'ChangeWTail'))==1
            
            %tails
            %fitP 4 to ... are priors
            for i = 1 : numPriors
                TailPrior(i) = fitP(3+i);
            end
            
            %k priors
            kpriors = fitP(11);
        end
    end
    
    %(case fat tail priors and llhs)
    if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))==1        
        if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==0
            
            %case tail fixed across priors, Kvm free
            if sum(strcmp(varargin{:},'ChangeWTail'))==0
                
                
                %tails
                TailPrior = fitP(11);
                
                %fitP 4 to ... are priors
                for i = 1 : numPriors
                    kpriors(i) = fitP(3+i);
                end
            end
            
            %case tail free across priors, Kvm fixed
            if sum(strcmp(varargin{:},'ChangeWTail'))==1
                
                %tails
                %fitP 4 to ... are priors
                for i = 1 : numPriors
                    TailPrior(i) = fitP(3+i);
                end
                
                %k priors
                kpriors = fitP(11);
            end
        end
    end
    
    %case fat tail prior and llh, all tails and Kvm free
    if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))==1
        if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==1
            
            %tails
            %WARNING: Make sure the fit parameters are in the following
            %order: kl24, kl12, kl6, kpriors1, kpriors2, kpriors3, kpriors4, kcard, Prand, Kmotor,
            %Tail_llh1, Tail_llh2, Tail_llh3, Tail_prior1, Tail_prior2, Tail_prior3, Tail_prior4
            
            %fitP numLLH to ... are priors
            for i = 1 : numPriors
                kpriors(i) = fitP(numLLH+i);
            end
                        
            %fitP 11 to ... are llh tails
            for i = 1 : numLLH
                                
                TailLLH(i) = fitP(numLLH + numPriors + 3 + i);
                
                %case Tail is fixed set TaiLLH = 0 if NaN
                if isnan(TailLLH(i))
                    TailLLH(i) = 0;
                end
            end
            
            %fitP 14 to ... are prior tails
            for i = 1 : numPriors
                TailPrior(i) = fitP(numLLH + numPriors + 3 + numLLH + i);
                
                %case Tail is fixed set TailPrior = 0 if NaN
                if isnan(TailPrior(i))
                    TailPrior(i) = 0;
                end
            end
        end
    end
    
    %(case von Mises without fat tail)
    if sum(strcmp(varargin{:},'FatTailPrior'))==0
        if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))==0
            
            %tails (should be NaN)
            TailPrior = NaN;
            
            %tails (should be NaN)
            TailLLH = zeros(3,1);
            
            %fitP 4 to ... are priors
            for i = 1 : numPriors
                kpriors(i) = fitP(numLLH + i);
            end
        end
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
end


%(case bimodal prior)
%--------------------
if sum(strcmp(priorShape,'bimodalPrior'))
    
    %modes and trials
    numPriors = size(SetOfpriorModes,2);
    for i = 1 : numPriors
        modesPrior(i,:) = SetOfpriorModes(i,:);
        Prior(:,i) = priorModes(:,1)==modesPrior(i,1) & priorModes(:,2)==modesPrior(i,2);
    end  
    
    %case prior can change tail's weight
    if sum(strcmp(varargin{:},'ChangeWTail'))==0
        
        %tails
        TailPrior = fitP(11);
        
        %fitP 4 to ... are priors
        for i = 1 : numPriors
            kpriors(i) = fitP(3+i);
        end
    end
    
    %case cannot
    if sum(strcmp(varargin{:},'ChangeWTail'))==1
        
        %tails
        %fitP 4 to ... are priors
        for i = 1 : numPriors
            TailPrior = fitP(3+i);
        end
        
        %k priors
        kpriors = fitP(11);
    end 
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

%motion directions
motdir = unique(feature);

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
        
        %coh 24,12,6
        [MAP,likeMAP1{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
        [~  ,likeMAP2{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
        [~  ,likeMAP3{i}] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
        
    end
elseif sum(strcmp(varargin{:},'SamplingReadout'))==1
    %(case Sampling readout)
    %case priors have same tail weight
    if numel(TailPrior)==1
        for i = 1 : numPriors
            %coh 24,12,6
            [MAP,likeMAP1{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl1,TailLLH(1),modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
            [~  ,likeMAP2{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl2,TailLLH(2),modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
            [~  ,likeMAP3{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl3,TailLLH(3),modesPrior(i,:),kpriors(i),kcardinal,TailPrior,priorShape,TheModel,varargin{:});
        end
    end
    
    %case priors have different tail weights
    if numel(TailPrior)>1
        for i = 1 : numPriors
            
            %coh 24,12,6
            [MAP,likeMAP1{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl1,TailLLH(1),modesPrior(i,:),kpriors(i),kcardinal,TailPrior(i),priorShape,TheModel,varargin{:});
            [~  ,likeMAP2{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl2,TailLLH(2),modesPrior(i,:),kpriors(i),kcardinal,TailPrior(i),priorShape,TheModel,varargin{:});
            [~  ,likeMAP3{i}] = SLBayesSamplingLookupTable(ulall,motdir,kl3,TailLLH(3),modesPrior(i,:),kpriors(i),kcardinal,TailPrior(i),priorShape,TheModel,varargin{:});
        end
    end
else
    fprintf('%s \n','(SLgetLoglBayesianModelcmae) The readout is missing. Please input either:')
    fprintf('%s \n','(SLgetLoglBayesianModelcmae) - "MAPReadout"')
    fprintf('%s \n','(SLgetLoglBayesianModelcmae) - "SamplingReadout"')
    keyboard
end

%--------
%%SET LLH
%--------
stimStrength24 = stimStrength==0.24;
stimStrength12 = stimStrength==0.12;
stimStrength06 = stimStrength==0.06;

%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials)
%for possible values of upo (rows) for each trial (column)
PupoGivenBI = nan(numel(MAP),numel(data));
for i = 1 : numel(motdir)
    
    %displayed motion direction for this trial
    thisd = feature==motdir(i);
    
    for j = 1 : numPriors
        
        %stimStrength 24,12,6
        PupoGivenBI(:,thisd&stimStrength24&Prior(:,j)) = likeMAP1{j}(:,motdir(i(ones(sum(thisd&stimStrength24&Prior(:,j)),1))));        
        PupoGivenBI(:,thisd&stimStrength12&Prior(:,j)) = likeMAP2{j}(:,motdir(i(ones(sum(thisd&stimStrength12&Prior(:,j)),1))));
        PupoGivenBI(:,thisd&stimStrength06&Prior(:,j)) = likeMAP3{j}(:,motdir(i(ones(sum(thisd&stimStrength06&Prior(:,j)),1))));
    end
end

%scale to probability.
%This is to make sure but it should already scale to probability and sum to
%1.
Z_ = sum(PupoGivenBI);
Z =  Z_(ones(size(PupoGivenBI,1),1),:);
PupoGivenBI =  PupoGivenBI./Z;

%probabilities of percepts "upo" given random estimation
PupoGivenRand =  ones(360,numel(feature))/360;

%calculate probability of percepts "upo" given the model
PBI =  1 - Prandom;
PupoGivenModel =  PupoGivenBI*PBI + PupoGivenRand*Prandom;

%check predicted densities
if any(isnan(sum(PupoGivenModel)))
    fprintf('%s \n','(SLgetLoglBayesianModelcmae) Predicted densities contained NaN.')
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    return
end

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel)) ==1
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
Pmotforconv = Pmot(:,ones(1,numel(feature)));
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
    fprintf('(SLgetLoglBayesianModelcmae) Some prediction densities were < but near 0, thus floored at 10^-320. \n')

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
    if sum(size(data)~=size([1:1:numel(feature)]'))==2
        data = data';
    end
    idx = sub2ind(size(PestimateGivenModel),data,[1:1:numel(feature)]');
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
        [meanData,stdData,dataCond] = SLmakeDataMeanAndStd(data,feature,...
            stimStrength,pstd,priorModes,priorShape);
        
        %pred
        [meanPred, stdPred] = makePred(feature,stimStrength,pstd,priorShape,priorModes,'Mean',[],PestimateGivenModel,MAP);
        
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
                
                fprintf('%.2f  %.2f   %.2f   %.2f   %.2f  %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.10f  %.2f \n',...
                    negLogl,kl1,kl2,kl3,kpriors,kcardinal,Prandom,km,TailPrior,ti)
            end
        end
            
            %case priors has heavy tail and same von mises strengths with different tail weight
        if sum(strcmp(varargin{:},'ChangeWTail'))==1
            
            fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.10f  %.2f  \n',...
                negLogl,kl1,kl2,kl3,TailPrior,kcardinal,Prandom,km,kpriors,ti)
        end
        
        %case priors and llh has free fat tails and von mises strengths
        if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==1
            
            fprintf('%s %.2f  %.2f %.2f %.2f  %.2f  %.2f  %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.10f  %.2f  \n',...
                '(SLgetLoglBayesianModel)',negLogl,kl1,kl2,kl3,kpriors,kcardinal,Prandom,km,TailLLH,TailPrior,ti)
        end
    end
end

%rearrange for cmaes
if strcmp(TheModel,'withoutCardinal')
   fitP(11) = [];
   fitP(8) = [];
end



%predictions
function [meanPred, stdPred] = makePred(feature,stimStrength,pstd,priorShape,priorModes,TrialOrMean,output,PestimateGivenModel,MAP)

%case von Mises prior
%--------------------
%conditions
if strcmp(priorShape,'vonMisesPrior')
    [cond,idxCondUniqtrial,~] = SLuniqpair([pstd stimStrength feature]);
end

%case bimodal prior
%------------------
%conditions
if strcmp(priorShape,'bimodalPrior')
    
    %check all prior conditions are there
    priorCond = priorModes(:,2) - priorModes(:,1);
    numPriorCond = size(SLuniqpair(priorModes),1);
    
    if numPriorCond == numel(unique(priorCond))
        
        [cond,idxCondUniqtrial,~] = SLuniqpair([priorCond stimStrength feature]);
        
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
            trialsThiC = pstd==cond(i,1) & stimStrength==cond(i,2) & feature==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
    
    %case bimodal prior
    %------------------
    if strcmp(priorShape,'bimodalPrior')
        output.TrialPred = nan(size(priorCond));
        for i = 1 : numCond
            trialsThiC = priorCond==cond(i,1) & stimStrength==cond(i,2) & feature==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
end