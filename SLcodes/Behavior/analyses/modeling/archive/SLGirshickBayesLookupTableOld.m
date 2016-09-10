% SLGirshickBayesLookupTable.m
%
%     author: steeve laquitaine
%       date: 140522 (last update 141121)
%
%    purpose: get the likelihood of circular data ranging from 1:1:360
%             given a Bayesian model with
%             
%             #a von Mises likelihood
%
%             #a learnt prior that can be 
%             - a von Mises v(modePrior,klearnt)
%               with modePrior the mode
%               k the concentration parameter
%
%             - a mixture of two von Mises v=0.5*v(modePrior(1),klearnt)+0.5*v(modePrior(2),klearnt)
%               with modePrior the two modes
%               k the concentration parameter
%
%             #a cardinal prior that is a mixture of 4 von mises with modes
%             at the cardinal directions and same concentration parameters
%             k
%
%     usage:
%
%       tic
%       [uniqallMAP,likelihoodUniqMAP] = SLGirshickBayesLookupTable(1:1:360,...
%       5:10:355,5,225,4.77,0,'vonMisesPrior','withCardinal','FatTailPrior');
%       toc
%
% Description: this is an implementation of the Girshick lookup table.
%              posterior = likelihood.*PRIORcardinal.*PRIORlearnt;
%              
%              The model accounts for cases when a same evidence could lead
%              two many equally possible percepts (>1 MAP percepts).
%               
%              duration: about 0.1 sec.    
%
%varargin
%--------
%
%
% Girshick et al., Nature Neuroscience.

function [uniqallMAP,likelihoodUniqMAP] = SLGirshickBayesLookupTableOld(ulall,...
    motdir,kl,modePrior,klearnt,kcardinal,priorShape,TheModel,varargin)

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
%likelihood of measurement mi(col) given motion directions di (row)
l = vmPdfs(diSpace,ulall,kl,'norm');


%---------------
%CARDNINAL PRIOR
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
if strcmp(varargin{:},'FatTailPrior')
    
end
    
%posteriors: probability of common causes explaining likelihood,
%cardinal and learnt priors. We fix probabilities at 10^10 floating
%points. This permits to the modes of the posterior despite round-off
%errors. Try with different combinations of 10^6 and round instead of fix
%If we don't round enough we cannot get the modes of the posterior
%accurately due to round-off errors. But, now if we
%round too much we get more modes than we should, but the values
%obtained surf around the values of the true modes so I choose to round
%more than not enough (same as in simulations).


%----------
%POSTERIORS
%----------
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
%find MAPs estimates (values) of each ul/mi (rows). whem a row  has two
%values it means that the same evidence produced two MAPs next to each other
%due to lack of numerical precision.
%e.g., when both likelihood and learnt prior are weak, an evidence produces
%a relatively flat posterior which maximum can not be estimated accurately.
%max(posterior) produces two (or more) consecutive MAP values. The best way
%is to record them and later average those MAP values to get the middle
%value.
MAP = nan(numel(m),length(diSpace));
for i = 1 : numel(m)
    numMAPs = sum(po(:,i)==max(po(:,i)));
    MAP(i,1:numMAPs) = diSpace(po(:,i)==max(po(:,i)));
end

%get the max likelihood of observing each data upo that is the probability
%of ul/mi given the motion direction for each ul and associated upo.
%ul/mi
maxnumMAP = max(sum(~isnan(MAP),2));
ul=m';
ul=ul(:,ones(maxnumMAP,1));
ul=ul(:);

%associate MAP estimate
MAP=MAP(:,1:maxnumMAP);
MAP=MAP(:);

%associate likelihood for each motion direction
likelihoodMAPgivenMi = repmat(mPdfs,maxnumMAP,1);

%everything sorted by MAP
TheMatrix = [ul MAP likelihoodMAPgivenMi];
TheMatrix = sortrows(TheMatrix,2);
ul = TheMatrix(:,1);
MAP = TheMatrix(:,2);
likelihoodMAPgivenMi = TheMatrix(:,3:end);


%----------------------------
%LIKELIHOOD EACH MAP ESTIMATE
%----------------------------

%old method (up to 140811)
%-------------------------
% %Set the likelihood of MAPs not produced at 0 because the model cannot
% %produce those estimates even at a reasonable resolutions of motion
% %direction.
% uniqMAP = unique(MAP(~isnan(MAP)));
% MAPnotProduced = setdiff(1:1:360,uniqMAP)';
% likelihoodMAPnotProduced = zeros(numel(MAPnotProduced),size(mPdfs,2));
% 
% %Add MAP not produced and their likelihood=0 & sort everything again by MAP
% likelihoodMAPs=[likelihoodMAPgivenMi;likelihoodMAPnotProduced];
% allMAP=[MAP;MAPnotProduced];
% ulnonExistent=nan(numel(MAPnotProduced),1);
% ul=[ul;ulnonExistent];
% TheMatrixII=[ul allMAP likelihoodMAPs];
% TheMatrixII=sortrows(TheMatrixII,2);
% allMAP=TheMatrixII(:,2);
% likelihoodMAPs=TheMatrixII(:,3:end);
% 
% %get the likelihood of each MAP (rows are MAPs, cols are motion directions,
% %values are likelihood). If a MAP had been produced by different
% %mi, then its likelihood is the average probability of this mi for each
% %motion direction, it is consistent with the laws of probability. It is
% %the (1/nummi)*P(MAP/mi1) + (1/nummi)*P(MAP/mi2) +...
% %It works very nice, the estimates density obtained are very smoothed.
% 
% %note: we can see horizontal stripes of 0 likelihood near the obliques when
% %cardinal prior is strong because obliques MAPs are never produced.
% %The range of MAPs not produced near the obliques increase significantly
% %with cardinal prior strength.
% uniqallMAP = unique(allMAP(~isnan(allMAP)));
% likelihoodUniqMAPtmp = nan(numel(uniqallMAP),size(mPdfs,2));
% for i = 1 : numel(uniqallMAP)
%     
%     %find mi that produced this same MAP
%     posUniqMAP = allMAP==uniqallMAP(i);
%     
%     %average probabilities over evidences mi that produces this same MAP
%     likelihoodUniqMAPtmp(i,:) = mean(likelihoodMAPs(posUniqMAP,:),1);
% end
% likelihoodUniqMAP = nan(numel(uniqallMAP),length(diSpace));
% likelihoodUniqMAP(:,motdir) = likelihoodUniqMAPtmp;
% 
% %to visualize
% %imagesc(likelihoodUniqMAPtmp)



%Best method (from 140811)
%-------------------------

%Set the likelihood of MAPs not produced at 0 because the model cannot
%produce those estimates even at a reasonable resolutions of motion
%direction.
uniqMAP = unique(MAP(~isnan(MAP)));
MAPnotProduced = setdiff(1:1:360,uniqMAP)';
numMissingMAP = numel(MAPnotProduced);
likelihoodMAPnotProduced = zeros(numMissingMAP,numMotdir);

%Add MAP not produced and their likelihood=0 & sort everything again by MAP
likelihoodMAPs = [likelihoodMAPgivenMi;likelihoodMAPnotProduced];
allMAP = [MAP;MAPnotProduced];
ulnonExistent = nan(numMissingMAP,1);
ul=[ul;ulnonExistent];
TheMatrixII = [ul allMAP likelihoodMAPs];
TheMatrixII = sortrows(TheMatrixII,2);
allMAP = TheMatrixII(:,2);
likelihoodMAPs = TheMatrixII(:,3:end);

%get the likelihood of each MAP (rows are MAPs, cols are motion directions,
%values are likelihood). If a MAP had been produced by different
%mi, then its likelihood is the average probability of this mi for each
%motion direction, it is consistent with the laws of probability. It is
%the (1/nummi)*P(MAP/mi1) + (1/nummi)*P(MAP/mi2) +...
%It works very nice, the estimates density obtained are very smoothed.

%note: we can see horizontal stripes of 0 likelihood near the obliques when
%cardinal prior is strong because obliques MAPs are never produced.
%The range of MAPs not produced near the obliques increase significantly
%with cardinal prior strength.
uniqallMAP = unique(allMAP(~isnan(allMAP)));
likelihoodUniqMAPtmp = nan(numel(uniqallMAP),numMotdir);
for i = 1 : numel(uniqallMAP)
    
    %find mi that produced this same MAP
    posUniqMAP = allMAP==uniqallMAP(i);
    
    %average probabilities over evidences mi that produces this same MAP
    likelihoodUniqMAPtmp(i,:) = nanmean(likelihoodMAPs(posUniqMAP,:),1);
end

%case some MAPs are not produced, Interpolate probability of observing 
%the missing MAPs
[~,PosObsMap] = setdiff(uniqallMAP,MAPnotProduced);

%There must be at least two MAPs produced to be able to interpolate
if ~isempty(MAPnotProduced) && numel(PosObsMap)>=2
    pmissing = interp1(uniqallMAP(PosObsMap),...
        likelihoodUniqMAPtmp(PosObsMap,:),MAPnotProduced,'linear');
    
    %now find MAPs that are most likely never produced. A MAP is most likely never
    %produced if there is a big gap between this missing MAP and the next 
    %closest MAP. It would be unlikely that this MAP would
    %be missing due to lack of numerical precision.
    %Then replace the likelihood of this missing MAP never produced by 0.
    for i = 1 : numMissingMAP
        NextClosestObsMAP = min(uniqMAP(uniqMAP > MAPnotProduced(i)));
        if ~isempty(NextClosestObsMAP)
            %This MAP is never produced if its distance to the next observed MAP is
            % > 3 deg. I tested different distance until I obtained smooth
            %distributions.
            if abs(vectors2signedAngle(MAPnotProduced(i),NextClosestObsMAP)) > 3;
                pmissing(i,:) = zeros(1,numMotdir);
            end
        else
            pmissing(i,:) = zeros(1,numMotdir);
        end
    end
    likelihoodUniqMAPtmp(MAPnotProduced,:) = pmissing;
    
elseif ~isempty(MAPnotProduced) && numel(PosObsMap)<2
    %case all MAPs are missing except one. Set the likelihood of all missing
    %MAP at zero (e.g., if the prior only is not flat).
    likelihoodUniqMAPtmp(MAPnotProduced,:) = zeros(numMissingMAP,numMotdir);
end
likelihoodUniqMAPtmp(isnan(likelihoodUniqMAPtmp)) = 0;

%place likelihood in a large matrix for all directions
likelihoodUniqMAP = nan(numel(uniqallMAP),length(diSpace));
likelihoodUniqMAP(:,motdir) = likelihoodUniqMAPtmp;

%to visualize
%imagesc(likelihoodUniqMAPtmp)