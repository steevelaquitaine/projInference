
% SLGirshickBayesLookupTable.m
%
%     author: steeve laquitaine
%       date: 140620 udpated 150131
%             last modification helped speed up the code a lot (>10x). We
%             only calculate the likelihood of the data given the motion
%             direction actually displayed (n=36) and not the full range 
%             1:1:360. logl result is identical to previous method but 10x
%             faster.
%    purpose: get the model likelihood of motion direction estimates.
%
%     usage:
%
%           [uniqallPoSmpl,likelihoodUniqPoSmpl] = SLBayesSamplingLookupTable(1:1:360,(5:10:355)',5,0.5,225,33,0,0.6,'vonMisesPrior','withoutCardinal','FatTailPrior');
%
%
% Description: A Bayesian sampling lookup table.
%              posterior = likelihood.*PRIORcardinal.*PRIORlearnt;
%              
%
% reference: Pouget et al.,PNAS, 2011

function [uniqallPoSmpl,likelihoodUniqPoSmpl] = SLBayesSamplingLookupTable(ulall,...
    motdir,kl,TailLLH,modePrior,klearnt,kcardinal,TailPrior,priorShape,TheModel,varargin)
%lookup table for Bayesian inference with cardinal and learnt prior


%call for help
% if ieNotDefined('ulall')
%     help SLBayesSamplingLookupTable
%     return
% end

%Sanity checks
if length(varargin{:})<2
    fprintf('%s \n','(SLBayesSamplingLookupTable) Somethings wrong. There',...
        'are not enough inputs in varargin...')
    keyboard
end

%To increase the resolution of PoSmpls estimates, we can increase the
%motion directions resolution "di" e.g. 1:0.5:360 instead of [1 1 360].
%It works fine.
%Direction space, displayed motion direction and evidences.
diSpace = (1:1:360)';
numMotdir = length(motdir);
m = 1:1:360;

%-------------------------
%MEASUREMENT DISTRIBUTIONS 
%-------------------------
%measurement probability densities ~v(diSpace; motion direction,kl) 
%di: displayed direction (col)
%over possible evidence m(row), typically 1:1:360.
mPdfs = vmPdfs(diSpace,motdir,kl,'norm');

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
    PRIORcardinal=vmPdfs(diSpace,[90 180 270 360],kcardinal,'norm');
    PRIORcardinal=0.25.*sum(PRIORcardinal,2);
    PRIORcardinal=PRIORcardinal(:,ones(numel(m),1));
end

%------------
%LEARNT PRIOR
%------------

%(case learnt prior is von Mises)
%------------------------------
%learnt prior (over motion directions (row) is same for each mi(col))
%mode is 225 degrees.
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
%modes are:
%[145 305]
%165 285]
%[185 265]
%[205 245]
if strcmp(priorShape,'bimodalPrior')
    
    %warning if two modes are not found
    if numel(modePrior)~=2
        fprintf('%s \n','(SLGirshickBayesLookupTable) There should be two modes not one..')
        keyboard
    end
    
    %learnt prior (over motion directions (row) is same for each mi(col))
    %repeat to combine with  each possible likelihood
    vm=vmPdfs(diSpace,modePrior,klearnt,'norm');
    PRIORlearnt=0.5*vm(:,1)+0.5*vm(:,2);
    PRIORlearnt=PRIORlearnt(:,ones(numel(m),1));
end

%(case learnt priors are fat tailed
%-----------------------------------------------------------------
if sum(strcmp(varargin{:},'FatTailPrior'))
        
    uniform = (1/360)*ones(size(PRIORlearnt));
    PRIORlearnt = (1-TailPrior)*PRIORlearnt + TailPrior*uniform;
    
end

%(case learnt prior & likelihood have same fat tail)
%---------------------------------------------------------
if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))

    %tail
    uniform = (1/360)*ones(size(PRIORlearnt));
    
    %prior and llh
    PRIORlearnt = (1 - TailPrior)*PRIORlearnt + TailPrior*uniform;
    l = (1 - TailPrior)*l + TailPrior*uniform;
end

%(case learnt prior & likelihood have different fat tails)
%---------------------------------------------------------
if sum(strcmp(varargin{:},'FatTailPriorAndLLH'))
    if sum(strcmp(varargin{:},'FitEachKvmAndTails'))
        
        %tail
        uniform = (1/360)*ones(size(PRIORlearnt));
        
        %fat tail prior
        PRIORlearnt = (1 - TailPrior)*PRIORlearnt + TailPrior*uniform;
        
        %fat tail LLH
        l = (1 - TailLLH)*l + TailLLH*uniform;
    end
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

%(case we fit cardinal priors)
%-----------------------------
%po (val=P(dir|evidence); row=dir, col=evidence)
if strcmp(TheModel,'withCardinal')
    po = l.*PRIORcardinal.*PRIORlearnt;
elseif strcmp(TheModel,'withoutCardinal')==1
    po = l.*PRIORlearnt;
end
Zpo = sum(po,1);
Zpo = Zpo(ones(numel(diSpace),1),:);
po = po./Zpo;
po = round(po*10^6)/10^6;

%(case von Mises prior)
%-----------------------
%TRICK: When k gets very large, e.g., for the prior, most values of the prior
%becomes 0 except close to the mean. The product of the likelihood and 
%prior only produces 0 values for all directions, particularly as motion
%direction is far from the prior. Marginalization (scaling) makes them NaN.
%If we don't correct for that fit is not possible. In reality a von Mises
%density will never have a zero value at its tails. We use the close-from
%equation derived in Murray and Morgenstern, 2010.
if strcmp(priorShape,'vonMisesPrior')
    
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


%Now get likelihood of data values (1:1:360) given the model
%find all estimates produced by each mi. A displayed motion direction
%can produce 360 different evidences. Each evidence can produce 360
%direction percepts from sampling.
%PoSmpl (nb mi  x  nb possible samples)
%This is the matrix that we want to create:
%
%   percepts     evidence      P(e|direction 1)  P(e|direction 2)   
%
%   sample 1  <-  ul = 1   <-  mpdf(1,dir=1)    mpdf(1,dir=2)   po
%   sample 1  <-     = 2   <-  mpdf(2,dir=1)    mpdf(1,dir=2)
%   ....      <-   
%   ....      <-     
%   sample 1  <-     = 360 <-  mpdf(360,dir=1)  mpdf(1,dir=2)
%
%   sample 2  <-  ul = 1   <-  mpdf(1,dir=1)        ...
%   sample 2  <-     = 2   <-  mpdf(2,dir=1)        ...
%   ....      <-  
%   ....      <-  .
%   sample 2  <-     = 360 <-  mpdf(360,dir=1)


%possibilities are 36 directions (5:10:355) x 360 measurements * 360 posterior samples
%   - P(measurement|di) is given by vm ~ (diSpace; di,km)
%   - P(sample|measurement) is given by posteriors
%directions
%motdir=1:360;
motdir = motdir';
motdirAll = motdir(ones(360*360,1),:);
motdirAll = motdirAll(:);

%measurements
% meas=1:360;
% meas=meas(ones(360,1),:);
% meas=meas(:);
% meas=meas(:,ones(360,1));
% meas=meas(:);

%samples
%360 possible sample percepts from the posterior
uniqallPoSmpl = (1:360)';
numSpl = length(uniqallPoSmpl);

%probability of evidence given true direction P(mi|di)
%mPdfs (val=p(mi|di), col=dir; row=evidence)
%"PmiGivenDi" is a column vector organized as follows:
%
% direc  evid  percept
%
% dir 1   mi1   sp1
%               sp2
%               ..
%               sp360
%
% dir 1   mi2   sp1
%               sp2
%               ..
%               sp360
%
% dir 2   mi1   sp1
%               sp2
%               ..
%               sp360

PmiGivenDi = SLreplicateRows(mPdfs,numSpl);
PmiGivenDi = PmiGivenDi(:);

%samples (36 di * 360 mi *360 sp)
%percept value of each row of "PmiGivenDi"
PoSmpl = uniqallPoSmpl(:,ones(length(motdir)*length(m),1));
PoSmpl = PoSmpl(:);

%p(sample|mi) i.e., posterior
%PsmpGivenMi = repmat(po(:),360,1);
%Bayesian posterior associated with each evidence mi.
%po (val=P(dir sample|evidence); row=sample, col=evidence)
potmp = po(:);
PsmpGivenMi = potmp(:,ones(numMotdir,1));
PsmpGivenMi = PsmpGivenMi(:);

%P(sample|dir,mi)
%PsampGivenDirAndMi = PmiGivenDi.*PsmpGivenMi;
PsampGivenDirAndMi = PmiGivenDi.*PsmpGivenMi;

%max llh(sample|dir) 
%It's: 
%P(sample=1|dir) =  P(sample=1|mi=1,dir)*P(mi=1|dir) + P(sample=1|mi=2,dir)*P(mi=2|dir) +...
%
%We loop over each direction (col) and each sample (row) and sum 
%the associated PsampGivenDirAndMi. This takes 13 seconds. It is the 
%fastest I could find. The function creates an array likelihoodUniqPoSmpl
%by accumulating elements of the vector PsampGivenDirAndMi using the 
%subscripts in [PoSmpl motdir];
likelihoodUniqPoSmpl = accumarray([PoSmpl motdirAll],PsampGivenDirAndMi);
dirNotDisplayed = setdiff(1:360,motdir);
likelihoodUniqPoSmpl(:,dirNotDisplayed) = NaN;

%see what it looks like
% imagesc(likelihoodUniqPoSmpl(:,motdir))





