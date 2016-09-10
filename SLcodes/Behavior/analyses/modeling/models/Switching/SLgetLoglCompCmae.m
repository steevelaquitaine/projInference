
%SLgetLoglCompCmae.m
%
%
% author: steeve laquitaine
%purpose: calculate logl data given switching model
%
%run with SLfitCompetitionModel
%
%
%inputs:
%
%       fitP : at least 10 params


function [negLogl,fitP,Logl_pertrial] = SLgetLoglCompCmae(fitP,data,...
    d,...
    stimStrength,...
    pstd,...
    priorShape,...
    priorModes,...
    TheModel,...
    varargin)
    
%get fit parameters. They should be in this order.
%(von Mises prior)
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal','Prand','km'
%(bimodal prior)
%'coh24','coh12','coh6','[145 305]','[165 285]','[185 265]','[205 245]','kcardinal','Prand','km'
%time
ticLogL = tic;

%model parameters
kl1 = fitP(1);
kl2 = fitP(2);
kl3 = fitP(3);
klearnt1 = fitP(4);
klearnt2 = fitP(5);
klearnt3 = fitP(6);
klearnt4 = fitP(7);

%no cardinal
if strcmp(TheModel,'withoutCardinal')
    kcardinal = 0;
    Prandom   = fitP(8);
    km        = fitP(9);
    %cardinal
elseif strcmp(TheModel,'withCardinal')
    kcardinal = fitP(8);
    Prandom   = fitP(9);
    km        = fitP(10);
end
fitP(8:10) = [kcardinal Prandom km];
ulall = 1:1:360;

%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time.
%Constrain other fraction parameters between 0 and 1, and other parameters
%as>0.
%and case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom > 1
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    return
end
if any(fitP < 0)
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
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

%(case learnt prior is a von Mises)
%----------------------------------
%if the prior is a von Mises, the unique mode of the prior(s) are those
%found in the databank from (SLMakedatabank).
%In our unimodal experiment it should be [225].
if sum(strcmp(priorShape,'vonMisesPrior'))
    
    %get the modes (should be 225)
    modesPrior1 = SetOfpriorModes(1,:);
    modesPrior2 = SetOfpriorModes(1,:);
    modesPrior3 = SetOfpriorModes(1,:);
    modesPrior4 = SetOfpriorModes(1,:);
    
    %find the trials for each prior
    Prior1 = pstd==80;
    Prior2 = pstd==40;
    Prior3 = pstd==20;
    Prior4 = pstd==10;
    
    %get strength of learnt prior
    klearnt(pstd==80) = fitP(4);
    klearnt(pstd==40) = fitP(5);
    klearnt(pstd==20) = fitP(6);
    klearnt(pstd==10) = fitP(7);
end


%(case learnt prior is bimodal)
%-----------------------------
%get the 4 priors in this order:

if sum(strcmp(priorShape,'bimodalPrior'))
    priorShape = 'bimodalPrior';
    
    %get the modes
    %[145 305]
    %165 285]
    %[185 265]
    %[205 245]
    modesPrior1 = SetOfpriorModes(1,:);
    modesPrior2 = SetOfpriorModes(2,:);
    modesPrior3 = SetOfpriorModes(3,:);
    modesPrior4 = SetOfpriorModes(4,:);
    
    %find the trials for each prior
    Prior1 = priorModes(:,1)==modesPrior1(1,1) & priorModes(:,2)==modesPrior1(1,2);
    Prior2 = priorModes(:,1)==modesPrior2(1,1) & priorModes(:,2)==modesPrior2(1,2);
    Prior3 = priorModes(:,1)==modesPrior3(1,1) & priorModes(:,2)==modesPrior3(1,2);
    Prior4 = priorModes(:,1)==modesPrior4(1,1) & priorModes(:,2)==modesPrior4(1,2);
    
    %The two von Mises that composes the bimodal prior have the same std 20?.
    klearnt(pstd==20) = fitP(6);
end

%-----------
%SET THE LLH
%-----------
%get stimStrengths and strength of evidence
%mean
motdir = unique(d);

%strength
coh24=stimStrength==0.24;
coh12=stimStrength==0.12;
coh06=stimStrength==0.06;
kl(stimStrength==0.24)=fitP(1);
kl(stimStrength==0.12)=fitP(2);
kl(stimStrength==0.06)=fitP(3);
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
%get likelihood of data in condition with stimStrength 24
[~,likeMAP11] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP12] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP13] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP14] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 12
[~,likeMAP21] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP22] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP23] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP24] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 6
[~,likeMAP31] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP32] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP33] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[MAP,likeMAP34] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

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

%--------------------------------------------------------
%Bayesian inference with learnt prior and cardinal priors
%--------------------------------------------------------
[~,likeMAP11] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP12] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP13] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP14] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 12
[~,likeMAP21] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP22] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP23] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP24] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 6
[~,likeMAP31] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP32] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP33] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[MAP,likeMAP34] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

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
weightPriorCardmn   = klearnt./(klearnt+kl);
weightLlhCardmn     = kl./(klearnt+kl);

%scale the mixing weights to probabilities (all sum to 1)
%this is not equal to 1. We want to make it equal to 1.
sumP=weightPriorCardmn + weightLlhCardmn + Prandom;
PpriorCardmn     = 1 - (Prandom + weightLlhCardmn)./sumP;
PllhCardmn       = 1 - (Prandom + weightPriorCardmn)./sumP;
Prandnew     = unique(1 - (weightLlhCardmn + weightPriorCardmn)./sumP);

%take just one value of Prandnew because sometimes because of numerical
%instability there is slight variation in the value of Prandnew which
%produces more than one value. Mathematically there should be only one
%value because priors and llh strength have been scaled some that they sum
%to 1 and Prandom is a constant. Thus Prandnew must be one value.
Prandnew=Prandnew(1);

%repeat P(choose prior), matrix values, of each trial (columns) for each
%possible data value (rows).
numm = numel(m);
PpriorCardmnall=PpriorCardmn(ones(numm,1),:);
PllhCardmnall=PllhCardmn(ones(numm,1),:);

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

%random choice and prior and like both have zero strength
PupoGivenModel(:,kl==0 & klearnt==0)=1/360;

%check predicted densities
if any(isnan(sum(PupoGivenModel)))
    fprintf('%s \n','(SLgetLoglCompCmae) Predicted densities contain NaN.')
    keyboard
end

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

%get read of extremely small negative numbers that occurs because of
%instability numerical instability. Set them to zero.
%PestimateGivenModel(PestimateGivenModel<0)=0;

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
upo=1:1:360;

%trial-predictions (sample estimate density). Variability in estimate
%should reflect variability in measurement density and motor noise
for i=1:length(d)
    pred(i)=randsample(1:1:360,1,'true',PestimateGivenModel(:,i));
end

%draw sampled predictions
%------------------------
%draw predictions and data
%drawCircStat(pred,d,coh,pstd);
%drawCircStat(data,d,coh,pstd);

%get loglikelihood of data
%-------------------------
%single trial's measurement, its position(row) for each trial(col) and its
%probability (also maxlikelihood of trial's data). Checked many times. It
%works.
%make sure sub2ind inputs are the same size
if sum(size(data)~=size([1:1:numel(d)]'))==2
    data=data';
end
idx = sub2ind(size(PestimateGivenModel),data,[1:1:numel(d)]');
PdataGivenModel = PestimateGivenModel(idx);

%We use log likelihood because likelihood is so small that matlab cannot
%encode it properly (numerical unstability). We can use single trials log
%likelihood to calculate AIC in the conditions that maximize differences in
%predictions of two models.
Logl_pertrial = log(PdataGivenModel);

%We use -sum(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
%logL = -sum(logL_pertrial);

%We use -sum(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
%negSumlogL=-sum(log(PdataGivenModel));

%We use E(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
negLogl = -sum(Logl_pertrial);

%Look at fitting. It is 3X faster without drawing.
ti = toc(ticLogL);
if sum(strcmp(varargin{1},'vonMisesPrior')) && sum(strcmp(varargin{1},'bimodalPrior'))
else
    fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
        negLogl,kl1,kl2,kl3,klearnt1,klearnt2,klearnt3,klearnt4,kcardinal,Prandom,km,ti)
end
