
% SLGirshickBayesLookupTable.m
%
%     author: steeve laquitaine
%       date: 140522 (last update 141121)
%
%    purpose: get the likelihood of circular data ranging from 1:1:360
%             given an augmented Bayesian model that predicts estimates (estimator) variability.
%
%      usage:
%
%               tic
%               [uniqallMAP,likelihoodUniqMAP] = SLGirshickBayesLookupTable(1:1:360,...
%               5:10:355,5,225,4.77,0,0,'vonMisesPrior','withCardinal','FatTailPrior');
%               toc
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
%              - a mixture of two von Mises v=0.5*v(modePrior(1),klearnt)+0.5*v(modePrior(2),klearnt)
%               with modePrior the two modes
%               k the concentration parameter
%
%               a cardinal prior that is a mixture of 4 von mises with modes
%               at the cardinal directions and same concentration parameters
%               k
%
%varargin
%--------
%
%
% Girshick et al., Nature Neuroscience.

function [uniqallMAP,likelihoodUniqMAP] = SLGirshickBayesLookupTable(ulall,...
    motdir,...
    kl,...
    modePrior,...
    klearnt,...
    kcardinal,...
    weightTail,...
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
diSpace = (1:1:360)';
numMotdir = length(motdir);
m = 1:1:360;

%-------------------------
%MEASUREMENT DISTRIBUTIONS 
%-------------------------
%measurement probability densities ~v(di,kl) 
%di the displayed direction (col)
%over range of possible measurement m(row), typically 1:1:360.
mPdfs = vmPdfs(m,motdir,kl,'norm');


%-----------
%LIKELIHOODS
%-----------
%likelihood of motion directions di (row) given measurements mi (col)
l = vmPdfs(diSpace,ulall,kl,'norm');


%---------------
%CARDINAL PRIOR
%---------------

%(case we fit cardinal priors)
%----------------------------
%cardinal prior (over motion directions (row) is same for each mi(col))
%case we don't want to fit the cardinal prior.
if strcmp(TheModel,'withCardinal')  && ~isnan(kcardinal)
    PRIORcardinal = vmPdfs(diSpace,[90 180 270 360],kcardinal,'norm');
    PRIORcardinal = 0.25.*sum(PRIORcardinal,2);
    PRIORcardinal = PRIORcardinal(:,ones(numel(m),1));
end


%------------
%LEARNT PRIOR
%------------

%(case learnt prior is von Mises)
%--------------------------------
%learnt prior (over motion directions (row) is same for each mi(col))
if strcmp(priorShape,'vonMisesPrior')
    PRIORlearnt = vmPdfs(diSpace,modePrior,klearnt,'norm');
    PRIORlearnt = PRIORlearnt(:,ones(numel(m),1));
end

%(case learnt prior is bimodal (a mixture of von Mises))
%-------------------------------------------------------
%modePrior contains the values of the two modes and here we assume that
%the strength of the two von Mises ccomposing the bimodal prior
%is the same (we could also assume that they are different) klearnt is
%a scalar. We also assume that the mixture coefficient is 0.5;
if strcmp(priorShape,'bimodalPrior')
    
    %warning if two modes are not found
    if numel(modePrior)~=2
        fprintf('%s \n',['(SLGirshickBayesLookupTable) There should be',...
        ' two modes not one..'])
        keyboard
    end
    
    %learnt prior (over motion directions (row) is same for each mi(col))
    %repeat to combine with  each possible likelihood
    vm = vmPdfs(diSpace,modePrior,klearnt,'norm');
    PRIORlearnt = 0.5*vm(:,1)+0.5*vm(:,2);
    PRIORlearnt = PRIORlearnt(:,ones(numel(m),1));
end

%(case learnt prior are fat tailed)
%----------------------------------
if sum(strcmp(varargin{:},'FatTailPrior'))

    uniform = (1/360)*ones(size(PRIORlearnt));
    PRIORlearnt = (1-weightTail)*PRIORlearnt + weightTail*uniform;
    
end
    
%(case learnt prior and likelihood are similarly fat tailed)
%-----------------------------------------------------------
if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))

    uniform = (1/360)*ones(size(PRIORlearnt));
    
    %prior and llh
    PRIORlearnt = (1 - weightTail)*PRIORlearnt + weightTail*uniform;
    l = (1 - weightTail)*l + weightTail*uniform;
end
                                        
%----------
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


%-------------
%MAP ESTIMATES
%-------------
%find MAPs estimates (values) of each mi (rows). 
%A mi may produce more than one MAP (more than one value for a row).
%e.g., when both likelihood and learnt prior are weak, an evidence produces
%a relatively flat posterior which maximum can not be estimated accurately.
%max(posterior) produces two (or more) MAP values. Those MAP values are
%observed with the same probability (e.g., 50% if two MAPs) given a mi.
MAP = nan(numel(m),length(diSpace));
for i = 1 : numel(m)
    numMAPs = sum(po(:,i)==max(po(:,i)));
    MAP(i,1:numMAPs) = diSpace(po(:,i)==max(po(:,i)));
end

%Likelihood of observing each data upo that is the P(mi|motion dir)
%of the mi that produced upo.
maxnumMAP = max(sum(~isnan(MAP),2));
ul = m';
ul = ul(:,ones(maxnumMAP,1));
ul = ul(:);

%associate MAP estimate
MAP = MAP(:,1:maxnumMAP);

%associate P(MAP|mi)
%[141122] - P(MAP|mi) is important because e.g., when the same mi produces a 
%posterior with two modes given a motion direction (flat likelihood and bimodal prior),
%the two modes are observed equiprobably P(MAP_1|mi) = P(MAP_2|mi) = 0.5 and P(MAP_1|di) = 
%P(MAP_1|di) = P(MAP_1|mi)*P(mi|di) = 0.5*P(mi|di) and in the same way 
%P(MAP_2|di) = P(MAP_2|mi)*P(mi|di) = 0.5*P(mi|di)
%The P(MAP|mi) of each MAP observed for a mi is 1/nb_of_MAP
numMAPeachmi = sum(~isnan(MAP),2);

%MAPs value only depends on mi value which determine likelihood position
%and thus posterior position. The displayed motion directions have not control on
%MAPs values. They only determine the Probability of observing MAPs.
%So the MAPs and the number of MAPs observed for each mi (row) is the same for 
%all motion directions displayed (col). Row (mi) repetitions appears in the matrix 
%when a mi produced a posterior with two modes.
%e.g., the matrices for a max number of MAP per mi of 2 
%
%   mPdfs_MAP1 . . dir_D
%       .
%       .
%      mi_M
%   mPdfs_MAP2 . . dir_D
%       .
%       .
%      mi_M
PMAPgivenMi = SLmakeColumn(1./numMAPeachmi);
PMAPgivenMi = PMAPgivenMi(:,ones(1,maxnumMAP));
PMAPgivenMi = PMAPgivenMi(:);
PMAPgivenMi = PMAPgivenMi(:,ones(1,numMotdir));

%associate P(mi|di) of each mi (row) for each motion dir di(col)
%e.g., the matrices for a max number of MAP per mi of 2 
%
%   mPdfs_MAP1 . . dir_D
%       .
%       .
%      mi_M
%   mPdfs_MAP2 . . dir_D
%       .
%       .
%      mi_M
PmiGivenDi = repmat(mPdfs,maxnumMAP,1);

%sort by MAP
TheMatrix = [ul MAP(:) PMAPgivenMi PmiGivenDi];
TheMatrix = sortrows(TheMatrix,2);
ul  = TheMatrix(:,1);
MAP = TheMatrix(:,2);
PMAPgivenMi = TheMatrix(:,3:2+numMotdir);
PmiGivenDi = TheMatrix(:,3+numMotdir:end);


%----------------------------
%LIKELIHOOD EACH MAP ESTIMATE
%----------------------------
%P(MAP|di) = P(MAP|mi)*P(mi|di)
PMAPgivenDi = PMAPgivenMi.* PmiGivenDi;

%until 141122
%PMAPgivenDi = PmiGivenDi;


%Best method (from 140811)
%-------------------------
%Set the likelihood of MAPs not produced at 0 because the model cannot
%produce those estimates even at a reasonable resolutions of motion
%direction.
uniqMAP = unique(MAP(~isnan(MAP)));
missingMAP = setdiff(1:1:360,uniqMAP)';
numMissingMAP = numel(missingMAP);
PmissingMAPgivenDi = zeros(numMissingMAP,numMotdir);

%Add MAP not produced and their likelihood=0 & sort everything again by MAP
PallMAPgivenDi = [PMAPgivenDi ; PmissingMAPgivenDi];
allMAP = [MAP ; missingMAP];

ulnonExistent = nan(numMissingMAP,1);
ul = [ul ; ulnonExistent];

TheMatrixII = [ul allMAP PallMAPgivenDi];
TheMatrixII = sortrows(TheMatrixII,2);
allMAP = TheMatrixII(:,2);
PallMAPgivenDi = TheMatrixII(:,3:end);

%likelihood of each MAP (rows are MAPs, cols are motion directions,
%values are likelihood). When a same MAP has been produced by different
%mi produced by the same motion direction, then the MAP's likelihood is 
%its mi likelihood.  The likelihoods are properly scaled at the end to sum
%to 1. e.g.,
%................................................................................
%   if only mi_1 produces?MAP=100º and mi_2 also produces?MAP=100º
% ?and mi_1 and mi_2 are both produced by the same motion direction di
%   P(MAP|di) = P(MAP|mi_1)*P(mi_1|di) + P(MAP|mi_2)*P(mi_2|dir)
%   P(MAP|mi_1) = P(MAP|mi_1) = 1 because both mi only produce one MAP (the same)
%   so P(MAP|di) = P(mi_1|di) + P(mi_2|dir)
%................................................................................
%
%note: we can see horizontal stripes of 0 likelihood near the obliques when
%cardinal prior is strong because obliques MAPs are never produced.
%The range of MAPs not produced near the obliques increase significantly
%with cardinal prior strength.
uniqallMAP = unique(allMAP(~isnan(allMAP)));
PuniqMAPgivenDi = nan(numel(uniqallMAP),numMotdir);

%find mi that produced this same MAP
%and average probabilities over evidences mi that produces this same MAP
for i = 1 : numel(uniqallMAP)
    
    posUniqMAP = allMAP==uniqallMAP(i);
    PuniqMAPgivenDi(i,:) = sum(PallMAPgivenDi(posUniqMAP,:),1);
    
    %until 141122
    %PuniqMAPgivenDi(i,:) = nanmean(PallMAPgivenDi(posUniqMAP,:),1);

end

% %% 
% %Interpolate the probability of observing missing MAPs 
% %(MAPs not produced)
% %There must be at least two MAPs produced to be able to interpolate
% [~,PosObsMap] = setdiff(uniqallMAP,missingMAP);
% if ~isempty(missingMAP) && numel(PosObsMap)>=2
%     pmissing = interp1(uniqallMAP(PosObsMap),...
%         PuniqMAPgivenDi(PosObsMap,:),missingMAP,'linear');
%     
%     %now find MAPs that are most likely never produced. A MAP is most likely never
%     %produced if there is a big gap between this missing MAP and the next 
%     %closest MAP. It would be unlikely that this MAP would
%     %be missing due to lack of numerical precision.
%     %Then replace the likelihood of this missing MAP never produced by 0.
%     for i = 1 : numMissingMAP
%         NextClosestObsMAP = min(uniqMAP(uniqMAP > missingMAP(i)));
%         if ~isempty(NextClosestObsMAP)
%             %This MAP is never produced if its distance to the next observed MAP is
%             % > 3 deg. I tested different distance until I obtained smooth
%             %distributions.
%             if abs(vectors2signedAngle(missingMAP(i),NextClosestObsMAP)) > 3;
%                 pmissing(i,:) = zeros(1,numMotdir);
%             end
%         else
%             pmissing(i,:) = zeros(1,numMotdir);
%         end
%     end
%     PuniqMAPgivenDi(missingMAP,:) = pmissing;
%     
% elseif ~isempty(missingMAP) && numel(PosObsMap)<2
%     %case all MAPs are missing except one. Set the likelihood of all missing
%     %MAP at zero (e.g., if the prior only is not flat).
%     PuniqMAPgivenDi(missingMAP,:) = zeros(numMissingMAP,numMotdir);
% end

PuniqMAPgivenDi(isnan(PuniqMAPgivenDi)) = 0;

%place likelihood in a large matrix for all 360 possible directions
likelihoodUniqMAP = nan(numel(uniqallMAP),length(diSpace));
likelihoodUniqMAP(:,motdir) = PuniqMAPgivenDi;

%scale to probability.
Z_ = sum(likelihoodUniqMAP);
Z =  Z_(ones(size(likelihoodUniqMAP,1),1),:);
likelihoodUniqMAP =  likelihoodUniqMAP./Z;







