
%SLsimSwitchingEstimatesDensity.m
%
% author: steeve laquitaine
%   date: 150807
%purpose: make Bayesian model motion estimates predictions.
%
%  usage:
%
%       [estDensity,estimates] = SLsimSwitchingEstimatesDensity(15,10,10,1e-3,83,0,'vonMisesPrior',225,'withoutCardinal','MAPReadout')
%
%
% inputs:
%
%        priorShape : 'vonMisesPrior' or 'bimodalPrior'
%          TheModel : 'withCardinal'
%          varargin : 'MAPReadout' or 'SamplingReadout'
%
% notes: 
%
% - most of the cases produce unimodal densities   
%
% - But for weak and distant llh and prior (e.g., 0.66 vs 0.66) : estimates densities are
%   bimodal. This is due to densities circularity.
%   Given a motion direction, evidence call fall withing 180 deg CW or within 180 deg CCW to
%   the prior producing a different peak in each case.
%

function [estDensity,estimates] = SLsimSwitchingEstimatesDensity(d,Kllh,Kprior,Prandom,km,kcardinal,priorShape,...
    priorModes,...
    TheModel,...
    varargin)

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

ulall = 1:1:360;

%--------------
%SET THE PRIORS
%--------------
%
%(case cardinal prior)
%---------------------
%case we don't want to fit the cardinal prior.
if isnan(kcardinal)
    
    TheModel = 'withoutCardinal';
    kcardinal = 0;
    
else
    
    %case we don't want to fit the cardinal prior.
    TheModel ='withCardinal';
    
end

%-----------
%SET THE LLH
%-----------

%---------------------------------------------------------
%Bayesian inference (BI) with evidence and cardinal priors
%---------------------------------------------------------
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
mykprior = 0;
mywTail  = 0;
[estimatesEC,likeMAPec] = SLGirshickBayesLookupTable(ulall,d,Kllh,priorModes,mykprior,kcardinal,mywTail,priorShape,TheModel,varargin{:});
likeMAPec = likeMAPec(~isnan(likeMAPec));
likeMAPec = likeMAPec./sum(likeMAPec);

%BI with learnt & cardinal priors
mykllh = 0;
[estimatesLC,likeMAPlc] = SLGirshickBayesLookupTable(ulall,d,mykllh,priorModes,Kprior,kcardinal,mywTail,priorShape,TheModel,varargin{:});
likeMAPlc = likeMAPlc(~isnan(likeMAPlc));
likeMAPlc = likeMAPlc./sum(likeMAPlc);

%check
if sum(estimatesEC ~= estimatesLC)==1
    fprintf('%s \n','(SLsimSwitchingEstimatesDensity) Estimate densities for prior and evidence components are not defined on the same estimate space')
    keyboard
end
estimates = estimatesLC;

%get p(percepts | "random process")
PupoGivenRand = ones(360,1)/360;

%calculate the fraction of trials that is controlled by each process 
%We use a divisive normalization competition rule between LLH, 
%learnt prior & cardinal prior. Basically the stronger the 
%representation and the
%most often it will be chosen. We use the strength of the von Mises
%composing the cardinal prior and not some measure of the overall strength
%of the cardinal prior. But both are correlated. For fixed modes, if the
%strength of the von Mises (one value here) increases, so does the strength
%of the overall distribution. I think...
%The probability of random estimation is fixed.

%case von Mises Prior
if sum(strcmp(priorShape,'vonMisesPrior')) == 1   

    wPriorCardmn = Kprior./(Kprior + Kllh);
    wLlhCardmn   = Kllh./(Kprior + Kllh);
    
end

%scale weights to sum to 1
sumP = wPriorCardmn + wLlhCardmn + Prandom;
PpriorCardmn  = 1 - (Prandom + wLlhCardmn)./sumP;
PllhCardmn    = 1 - (Prandom + wPriorCardmn)./sumP;
Prandnew      = unique(1 - (wLlhCardmn + wPriorCardmn)./sumP);

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
PupoGivenModel = likeMAPec.*PllhCardmn+ likeMAPlc.*PpriorCardmn + PupoGivenRand.*Prandnew;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel))==1
    keyboard
    fprintf('%s \n','Something wrong. PupoGivenModel are proba. thus must sum to 1')
end

%convolve with motor noise
%--------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo = 0 : 1 : 359;
Pmot = vmPdfs(upo,0,km,'norm');
estDensity = SLcircConv(PupoGivenModel,Pmot);

%get read of extremely small negative numbers that occurs because of 
%instability numerical instability. Set them to zero.
if ~isempty(estDensity<=0)
    
    estDensity(estDensity<=0) = 10^-320;
    fprintf('(SLgetLoglBayesianModel) P(estimates|model) were sometimes <0 but near zero (< 10^-10),thus floored at the lowest positive numerical value 10^-320. \n')
    
end

%normalize to probability; If we don't normalize, more random
%choice always increase probability of observing data causing larger
%probability of random choice to prevail. To avoid that we need to normalize
%to probabilities that sum to 1. It also makes intuitive sense to deal
%with probabilities.
%Here too estDensity integrated over estimates should already sum to
%1. But we still scale to be absolutely sure.
estDensity = estDensity./sum(estDensity);