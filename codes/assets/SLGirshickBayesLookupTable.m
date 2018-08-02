
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
%               [uniqallpercept,likelihoodUniqpercept] = SLGirshickBayesLookupTable(1:1:360,...
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
%              two many equally possible percepts (>1 percept percepts).
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

function [uniqallpercept,likelihoodUniqpercept] = SLGirshickBayesLookupTable(ulall,...
    motdir,...
    kl,...
    modePrior,...
    klearnt,...
    kcardinal,...
    weightTail,...
    priorShape,...
    TheModel,...
    varargin)

%lookup table for Bayesian inference with cardinal and learnt prior
%To increase the resolution of percepts estimates, we can increase the
%motion directions resolution "di" e.g. 1:0.5:360 instead of [1 1 360].
%It works fine.
%We calculate likelihood of data only for the 36 directions actually
%displayed and not the full range of 360 directions. This speed up the
%code by 10 times. Direction space, displayed motion direction and evidences.
diSpace = (1:1:360)';
numMotdir = length(motdir);
m = 1:1:360;

%case motion energy model
if slIsInput(varargin{:},'motionenergymodel')
    motionStrength = varargin{:}{find(strcmp(varargin{:},'motionStrength'))+1};
    %load precomputed sensory likelihood matrix for the current motion
    %strength
    if motionStrength==0.06
        load likelihood006.mat
        l = likelihood006;
    elseif motionStrength==0.12
        load likelihood012.mat
        l = likelihood012;
    elseif motionStrength==0.24
        load likelihood024.mat
        l = likelihood024;
    end
    %evidence distribution matrix is the matrix transpose of the sensory
    %likelihood
    mPdfs = l(motdir,:)';
end

%case it is another model, compute evidence distributions and
%sensory likelihood with von Mises shapes.
if ~slIsInput(varargin{:},'motionenergymodel')
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
end

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
% if strcmp(priorShape,'bimodalPrior')
%
%     %warning if two modes are not found
%     if numel(modePrior)~=2
%         fprintf('%s \n',['(SLGirshickBayesLookupTable) There should be',...
%         ' two modes not one..'])
%         keyboard
%     end
%
%     %learnt prior (over motion directions (row) is same for each mi(col))
%     %repeat to combine with  each possible likelihood
%     vm = vmPdfs(diSpace,modePrior,klearnt,'norm');
%     PRIORlearnt = 0.5*vm(:,1)+0.5*vm(:,2);
%     PRIORlearnt = PRIORlearnt(:,ones(numel(m),1));
% end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           BAYESIAN INTEGRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% %find evidence (cols of Likelihood matrix "l")
% %that produces null posteriors
% %calculate the posterior mean (same as mode for
% %such strong posterior) and make those evidence
% %posteriors delta function (very strong posteriors)
% %that peak at those means.
% if ~all(Zpo>0)
%     m_po0 = m(~(Zpo>0));
%     mpr = SLde2r(modePrior,0);
%     me = SLde2r(m_po0,0);
%     upo_po0 = round(SLra2d(mpr + atan2( sin(me-mpr) , klearnt/kl+cos(me-mpr) )));
%     upo_po0(upo_po0==0)=360;
%     for i = 1 : sum(~(Zpo>0))
%         po(upo_po0(i),m_po0(i)) = 1;
%     end
% end
% Zpo = sum(po,1);
% if ~all(Zpo>0)
%     dbstack
%     fprintf('%s \n','(SLGirshickBayesLookupTable) all posteriors must sum to 1')
%     keyboard
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  READOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(case von Mises prior)
%-----------------------
%TRICK: When k gets very large, e.g., for the prior, most values of the prior
%becomes 0 except close to the mean. The product of the likelihood and
%prior only produces 0 values for all directions, particularly as motion
%direction is far from the prior. Marginalization (scaling) makes them NaN.
%If we don't correct for that fit is not possible. In reality a von Mises
%density will never have a zero value at its tails. We use the close-from
%equation derived in Murray and Morgenstern, 2010.
if ~slIsInput(varargin{:},'motionenergymodel') && strcmp(priorShape,'vonMisesPrior')
    
    %find the trials where it happens and correct
    pos=find(isnan(po(1,:)));
    if ~isempty(pos)
        %closed-form equation (Murray and Morgenstern., 2010)
        %mode of posterior
        mi=diSpace(pos);
        mirad=SLde2r(mi,1);
        uprad=SLde2r(modePrior,1);
        upo = round(mirad + atan2(sin(uprad - mirad),(kl./klearnt)' + cos(uprad - mirad)));
        %make sure upo belongs to diSpace
        upo = round(SLra2d(upo)');
        upo(upo==0) = 360;
        upo(upo>360) = 360; %should never happen
        upo(upo<1) = 360;%should never happen
        
        %strength of posterior
        kpo=sqrt(kl^2 + klearnt^2 + 2*klearnt.*kl.*cos(uprad - mirad));
        
        %create those posterior
        po(:,pos)=vmPdfs(diSpace,upo,kpo,'norm');
    end
end

%-----------
%MAP READOUT
%-----------
%find MAPs estimates (values) of each mi (rows).
%A mi may produce more than one MAP (more than one value for a row).
%e.g., when both likelihood and learnt prior are weak, an evidence produces
%a relatively flat posterior which maximum can not be estimated accurately.
%max(posterior) produces two (or more) MAP values. Those MAP values are
%observed with the same probability (e.g., 50% if two MAPs) given a mi.
if slIsInput(varargin{:},'MAPReadout')
    percept = nan(numel(diSpace),length(m));
    for i = 1 : numel(m)
        MAPpos = po(:,i)==max(po(:,i));
        numMAPs = sum(MAPpos);
        percept(i,1:numMAPs) = diSpace(MAPpos);
    end
    %-----------
    %BLS READOUT
    %-----------
elseif slIsInput(varargin{:},'BLSReadout')
    percept = nan(numel(diSpace),length(m));
    for i = 1 : length(m)
        %percept(i,1) = round(sum(po(:,i).*diSpace));
        tmp = SLcircWeightedMeanStd(diSpace,po(:,i));
        percept(i,1) = round(tmp.deg.mean);
    end
else
    fprintf('(SLGirshickBayesLookupTable) Readout is missing - Even Switching must be MAPReadout or BLSReadout')
    keyboard        
end
percept(percept==0) = 360;


%sanity check
if any(percept(:) > 360) || any(percept(:) < 1)
    fprintf('(SLGirshickBayesLookupTable) percept must be >=1 and <=360')
    keyboard
end


%Likelihood of observing each data upo that is the P(mi|motion dir)
%of the mi that produced upo.
maxnumpercept = max(sum(~isnan(percept),2));
ul = m';
ul = ul(:,ones(maxnumpercept,1));
ul = ul(:);

%associate percept estimate
percept = percept(:,1:maxnumpercept);

%associate P(percept|mi)
%[141122] - P(percept|mi) is important because e.g., when the same mi produces a
%posterior with two modes given a motion direction (flat likelihood and bimodal prior),
%the two modes are observed equiprobably P(percept_1|mi) = P(percept_2|mi) = 0.5 and P(percept_1|di) =
%P(percept_1|di) = P(percept_1|mi)*P(mi|di) = 0.5*P(mi|di) and in the same way
%P(percept_2|di) = P(percept_2|mi)*P(mi|di) = 0.5*P(mi|di)
%The P(percept|mi) of each percept observed for a mi is 1/nb_of_percept
numpercepteachmi = sum(~isnan(percept),2);

%percepts value only depends on mi value which determine likelihood position
%and thus posterior position. The displayed motion directions have not control on
%percepts values. They only determine the Probability of observing percepts.
%So the percepts and the number of percepts observed for each mi (row) is the same for
%all motion directions displayed (col). Row (mi) repetitions appears in the matrix
%when a mi produced a posterior with two modes.
%e.g., the matrices for a max number of percept per mi of 2
%
%   mPdfs_percept1 . . dir_D
%       .
%       .
%      mi_M
%   mPdfs_percept2 . . dir_D
%       .
%       .
%      mi_M
PperceptgivenMi = SLmakeColumn(1./numpercepteachmi);
PperceptgivenMi = PperceptgivenMi(:,ones(1,maxnumpercept));
PperceptgivenMi = PperceptgivenMi(:);
PperceptgivenMi = PperceptgivenMi(:,ones(1,numMotdir));

%associate P(mi|di) of each mi (row) for each motion dir di(col)
%e.g., the matrices for a max number of percept per mi of 2
%
%   mPdfs_percept1 . . dir_D
%       .
%       .
%      mi_M
%   mPdfs_percept2 . . dir_D
%       .
%       .
%      mi_M
PmiGivenDi = repmat(mPdfs,maxnumpercept,1);

%sort by percept
TheMatrix = [ul percept(:) PperceptgivenMi PmiGivenDi];
TheMatrix = sortrows(TheMatrix,2);
ul  = TheMatrix(:,1);
percept = TheMatrix(:,2);
PperceptgivenMi = TheMatrix(:,3:2+numMotdir);
PmiGivenDi = TheMatrix(:,3+numMotdir:end);


%----------------------------
%LIKELIHOOD EACH percept ESTIMATE
%----------------------------
%P(percept|di) = P(percept|mi)*P(mi|di)
PperceptgivenDi = PperceptgivenMi.* PmiGivenDi;

%until 141122
%PperceptgivenDi = PmiGivenDi;


%Best method (from 140811)
%-------------------------
%Set the likelihood of percepts not produced at 0 because the model cannot
%produce those estimates even at a reasonable resolutions of motion
%direction.
uniqpercept = unique(percept(~isnan(percept)));
missingpercept = setdiff(1:1:360,uniqpercept)';
numMissingpercept = numel(missingpercept);
PmissingperceptgivenDi = zeros(numMissingpercept,numMotdir);

%Add percept not produced and their likelihood=0 & sort everything again by percept
PallperceptgivenDi = [PperceptgivenDi ; PmissingperceptgivenDi];
allpercept = [percept ; missingpercept];

ulnonExistent = nan(numMissingpercept,1);
ul = [ul ; ulnonExistent];

TheMatrixII = [ul allpercept PallperceptgivenDi];
TheMatrixII = sortrows(TheMatrixII,2);
allpercept = TheMatrixII(:,2);
PallperceptgivenDi = TheMatrixII(:,3:end);

%likelihood of each percept (rows are percepts, cols are motion directions,
%values are likelihood). When a same percept has been produced by different
%mi produced by the same motion direction, then the percept's likelihood is
%its mi likelihood.  The likelihoods are properly scaled at the end to sum
%to 1. e.g.,
%................................................................................
%   if only mi_1 produces?percept=100? and mi_2 also produces?percept=100?
% ?and mi_1 and mi_2 are both produced by the same motion direction di
%   P(percept|di) = P(percept|mi_1)*P(mi_1|di) + P(percept|mi_2)*P(mi_2|dir)
%   P(percept|mi_1) = P(percept|mi_1) = 1 because both mi only produce one percept (the same)
%   so P(percept|di) = P(mi_1|di) + P(mi_2|dir)
%................................................................................
%
%note: we can see horizontal stripes of 0 likelihood near the obliques when
%cardinal prior is strong because obliques percepts are never produced.
%The range of percepts not produced near the obliques increase significantly
%with cardinal prior strength.
uniqallpercept = unique(allpercept(~isnan(allpercept)));
PuniqperceptgivenDi = nan(numel(uniqallpercept),numMotdir);

%find mi that produced this same percept
%and average probabilities over evidences mi that produces this same percept
for i = 1 : numel(uniqallpercept)
    
    posUniqpercept = allpercept==uniqallpercept(i);
    PuniqperceptgivenDi(i,:) = sum(PallperceptgivenDi(posUniqpercept,:),1);
    
    %until 141122
    %PuniqperceptgivenDi(i,:) = nanmean(PallperceptgivenDi(posUniqpercept,:),1);
    
end

% %%
% %Interpolate the probability of observing missing percepts
% %(percepts not produced)
% %There must be at least two percepts produced to be able to interpolate
% [~,PosObspercept] = setdiff(uniqallpercept,missingpercept);
% if ~isempty(missingpercept) && numel(PosObspercept)>=2
%     pmissing = interp1(uniqallpercept(PosObspercept),...
%         PuniqperceptgivenDi(PosObspercept,:),missingpercept,'linear');
%
%     %now find percepts that are most likely never produced. A percept is most likely never
%     %produced if there is a big gap between this missing percept and the next
%     %closest percept. It would be unlikely that this percept would
%     %be missing due to lack of numerical precision.
%     %Then replace the likelihood of this missing percept never produced by 0.
%     for i = 1 : numMissingpercept
%         NextClosestObspercept = min(uniqpercept(uniqpercept > missingpercept(i)));
%         if ~isempty(NextClosestObspercept)
%             %This percept is never produced if its distance to the next observed percept is
%             % > 3 deg. I tested different distance until I obtained smooth
%             %distributions.
%             if abs(vectors2signedAngle(missingpercept(i),NextClosestObspercept)) > 3;
%                 pmissing(i,:) = zeros(1,numMotdir);
%             end
%         else
%             pmissing(i,:) = zeros(1,numMotdir);
%         end
%     end
%     PuniqperceptgivenDi(missingpercept,:) = pmissing;
%
% elseif ~isempty(missingpercept) && numel(PosObspercept)<2
%     %case all percepts are missing except one. Set the likelihood of all missing
%     %percept at zero (e.g., if the prior only is not flat).
%     PuniqperceptgivenDi(missingpercept,:) = zeros(numMissingpercept,numMotdir);
% end

PuniqperceptgivenDi(isnan(PuniqperceptgivenDi)) = 0;

%place likelihood in a large matrix for all 360 possible directions
likelihoodUniqpercept = nan(numel(uniqallpercept),length(diSpace));
likelihoodUniqpercept(:,motdir) = PuniqperceptgivenDi;

%scale to probability.
Z_ = sum(likelihoodUniqpercept);
Z =  Z_(ones(size(likelihoodUniqpercept,1),1),:);
likelihoodUniqpercept =  likelihoodUniqpercept./Z;







