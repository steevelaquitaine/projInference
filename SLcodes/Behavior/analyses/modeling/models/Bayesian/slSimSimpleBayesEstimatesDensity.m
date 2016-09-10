
%slSimSimpleBayesEstimatesDensity.m
%
% author: steeve laquitaine
%   date: 150805
%purpose: make Bayesian model motion estimates density predictions.
%
%  usage:
%
%       [estDensity,estimates] = slSimSimpleBayesEstimatesDensity(15,1.89,0.66,1e-3,83,0,'vonMisesPrior',225,'withoutCardinal','MAPReadout')
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

function [estDensity,estimates] = slSimSimpleBayesEstimatesDensity(d,Kllh,Kprior,Prandom,km,kcardinal,priorShape,...
    priorModes,...
    TheModel,...
    varargin)
 
ulall = 1 : 1 : 360;

%------
%PRIORS
%------
%cardinal priors
if isnan(kcardinal)
    kcardinal = 0;
end

%1.Bayesian inference
%get a matrix (360 possible MAPs, 360 possible motion directions) of
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.

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
%to Kprior3 in the von Mises prior case).

%In the case of bimodal prior, the fitting procedure uses only Klearnt3
%the 6th fit parameter to estimate the four prior strength (whatever what
%we input in the 3 other prior strength parameters).

%In the case of von Mises prior, the four Kprior1 to Kprior 3 parameters
%are used to fit the data.

%In the case of combined prior, the three Kprior1, Kprior2, Kprior4
%parameters are used to fit vonMises data, and Kprior3 is used to fit
%both for vonMises and bimodal prior data.
if sum(strcmp(varargin{:},'MAPReadout'))
    
    [estimates,likeMAP1] = SLGirshickBayesLookupTable(ulall,d,Kllh,priorModes,Kprior,kcardinal,0,priorShape,TheModel,varargin{:});
    
end

%(case Sampling readout)
%-----------------------
if sum(strcmp(varargin{:},'SamplingReadout'))==1
    
    [estimates,likeMAP1] = SLBayesSamplingLookupTable(ulall,d,Kllh,0,priorModes,Kprior,kcardinal,0,priorShape,TheModel,varargin{:});
    
end
likeMAP1 = likeMAP1(~isnan(likeMAP1));

%scale to probability.
likeMAP1 =  likeMAP1./sum(likeMAP1);
PestGivenRand = 1/360; 
percDensity = likeMAP1*(1 - Prandom) + PestGivenRand*Prandom; %add random estimation

%check PupoGivenModel sum to 1
if ~unique(sum(percDensity )) == 1
    keyboard
    fprintf('%s \n','Something wrong. "estDensity" are probabilities and should sum to 1')
end

%convolve with motor noise
%-------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo = 0 : 1 : 359;
Pmot = vmPdfs(upo,0,km,'norm');
estDensity = SLcircConv(percDensity,Pmot);

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

if ~isempty(find(estDensity<=0))
    
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

