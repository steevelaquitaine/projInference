% slsimulateBayesianModelCircularEstimateSpace.m
%
%     author: steeve laquitaine
%       date: 140522 (last update 141121)
%
%    purpose: get the likelihood of circular data ranging from 1:1:360
%             given an augmented Bayesian model that predicts estimates (estimator) variability.
%
%      usage:
%
%
%       [po,MAP] = slsimulateBayesianModelCircularEstimateSpace(5:20:360,5:10:360,5,225,4.77,0,'vonMisesPrior','withoutCardinal');
%
%
%
%input: 
%
%         motdir :
%             kl :
%      modePrior :
%        klearnt :
%      kcardinal : 
%     weightTail :  0: no tail to 1: flat
%     priorShape : 'vonMisesPrior' or 'bimodalPrior'
%       TheModel : 'withCardinal' or 'withoutCardinal'
%       varargin :
%
% 
% Description: this is an implementation of the Girshick lookup table.
%              posterior = likelihood.*PRIORcardinal.*PRIORlearnt;
%              
%              The model accounts for cases when a same evidence could lead
%              two many equally possible percepts (>1 MAP percepts). 
%
%              a von Mises likelihood
%
%              a learnt prior that can be 
%              - a von Mises v(modePrior,klearnt)
%               with modePrior the mode
%               k the concentration parameter
%
%
% Girshick et al., Nature Neuroscience.

function [po,MAP] = slsimulateBayesianModelCircularEstimateSpace(ulall,...
    motdir,...
    kl,...
    modePrior,...
    klearnt,...
    kcardinal,...
    priorShape,...
    TheModel,...
    varargin)

%Sanity checks
% if length(varargin{:})<1
%     fprintf('%s \n','(SLGirshickBayesLookupTable) Somethings wrong. There',...
%         'are not enough inputs in varargin...')
%     keyboard
% end

%lookup table for Bayesian inference with cardinal and learnt prior
%To increase the resolution of MAPs estimates, we can increase the
%motion directions resolution "di" e.g. 1:0.5:360 instead of [1 1 360].
%It works fine.
%We calculate likelihood of data only for the 36 directions actually 
%displayed and not the full range of 360 directions. This speed up the 
%code by 10 times. Direction space, displayed motion direction and evidences.
diSpace = (5:10:360)';
numMotdir = length(motdir);
m = 5:10:360;
m = ulall;

%MEASUREMENT DISTRIBUTIONS 
%measurement probability densities ~v(di,kl) 
%di the displayed direction (col)
%over range of possible measurement m(row), typically 1:1:360.
%mPdfs = vmPdfs(m,motdir,kl,'norm');
mPdfs = vmPdfs(ulall,motdir,kl,'norm');


%LIKELIHOODS
%likelihood of motion directions di (row) given measurements mi (col)
l = vmPdfs(diSpace,ulall,kl,'norm');


%CARDINAL PRIOR
%cardinal prior (over motion directions (row) is same for each mi(col))
%case we don't want to fit the cardinal prior.
if strcmp(TheModel,'withCardinal')  && ~isnan(kcardinal)
    PRIORcardinal = vmPdfs(diSpace,[90 180 270 360],kcardinal,'norm');
    PRIORcardinal = 0.25.*sum(PRIORcardinal,2);
    PRIORcardinal = PRIORcardinal(:,ones(numel(m),1));
end


%(case learnt prior is von Mises)
%learnt prior (over motion directions (row) is same for each mi(col))
if strcmp(priorShape,'vonMisesPrior')
    PRIORlearnt = vmPdfs(diSpace,modePrior,klearnt,'norm');
    PRIORlearnt = PRIORlearnt(:,ones(numel(m),1));
end
                                        
%POSTERIORS
%----------
%posteriors: probability of common causes explaining likelihood,
%cardinal and learnt priors. We fix probabilities at 10^10 floating
%points. This permits to get the modes of the posterior despite round-off
%errors. Try with different combinations of 10^6 and round instead of fix
%If we don't round enough we cannot get the modes of the posterior
%accurately due to round-off errors. But, now if we
%round too much we get more modes than we should, but the values
%obtained surf around the values of the true modes so I choose to round
%more than not enough (same as in simulations).
if strcmp(TheModel,'withCardinal')
    po = l.*PRIORcardinal.*PRIORlearnt;    
elseif strcmp(TheModel,'withoutCardinal')==1
    po = l.*PRIORlearnt;
end
Zpo = sum(po,1);
Zpo = Zpo(ones(numel(diSpace),1),:);
po = po./Zpo;
po = round(po*10^6)/10^6;

%Estimate space
MAP = nan(numel(m),length(diSpace));
for i = 1 : numel(m)
    numMAPs = sum(po(:,i)==max(po(:,i)));
    MAP(i,1:numMAPs) = diSpace(po(:,i)==max(po(:,i)));
end

%plot posterior and estimate space
colors = linspecer(length(ulall));
for i = 1 : length(ulall)
    hold all; plot(diSpace,po(:,i),'color',colors(i,:))
    plot(MAP(i,1),zeros(1,length(MAP)),'.','markersize',30,'color',colors(i,:))
end
xlim([0 360])