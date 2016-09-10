




%slMakePredictionsWJMlgTprior.m
%
%
% author: steeve laquitaine
%   date: 150724
%purpose: calculate estimates predictions with WJma
%         model (global likelihood integrates local dot likelihood)
%         and the logl of data given the model.
%
%  usage:
%
%
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJMlgTprior({'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},[4 500 500 500 500 0.1 60 0.4],'experiment','vonMisesPrior')
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJMlgTprior({'sub01','sub08'},[4 500 500 500 500 0.1 60 0.4],'experiment','vonMisesPrior')
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJMlgTprior({'sub01'},[4 500 500 500 500 0.1 60 0.4],'experiment','vonMisesPrior','loadData')
%
%
%list of analyses
%----------------
%            subjects : e.g., {'sub01'} or e.g., {'sub01','sub02'}
%           'loadData': load existing database
%
%
%
%Description:
%
%       - Motion direction likelihood is the integrated likelihood densities of the local dots directions
%       - Bayesian integration of Like and prior into posterior
%       - Estimates are posterior mean
%
%
%notes:
%   - 12 cores: duration: 25 sec (5000 sim trials)


function [Pdata,pPred,Bins,o] = slMakePredictionsWJMlgTprior(subjects,simp,varargin)

%Prealloc
global x
global coherent
global like
global su24
global su12
global su6
global Nsu24
global Nsu12
global Nsu6
global Ntrials
global Ndots
global Ns
global nCon
global svec
global pmean
global coh
global prior
global predEsttmp
global pos
global tailtmp
global amp1tmp
global rotrad
global motevec
global su
global Nsu
global sPest
global prPest
global cohPest
global w
global M
global ix
global Mcond
global rotdeg
global shifttmp

%--------
%get data
%--------
fprintf('%s \n','(slMakePredictionsWJMlgTprior) Gathering data ...')
fprintf('%s \n','(slMakePredictionsWJMlgTprior) Please set dataPath ...')

%load data from path
if any(strcmp(varargin,'loadData'))
    fprintf('(slMakePredictionsWJMlgTprior) Loading data from directory.')
    load datbank
    db = databank;
else
    %dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');
    dataPath = '~/data/dataPsychophy/proj01_priorStrength';
    %dataPath = '~/Dropbox/myDropbox/Project_withJustin/data/datapsychophy/proj01_priorstrength';
    cd(dataPath)
    varg = [varargin,'dataPath',dataPath];
    db   = SLMakedatabank(subjects,varg);
    fprintf('%s \n','Making database from',dataPath)
end

%gather data
data           = db.estimatesDeg;     %estimates
sStrg          = db.stimStrength;     %coherence
s              = db.stimFeatureDeg;   %direction
priorModes     = [db.priormodes{:}];  %prior modes
pstd           = db.Pstd;             %prior strength

%data densities by bin by 202 conditions
dataBins       = 0:10:360;
[Pdata,xpdf,o] = slmakeDataDistOverSub(subjects,db,dataBins,varargin{:});
o.posC         = o.posC;
nCon           = o.nCond;
coh            = sort(unique(o.uniqCond(:,2)),'descend');
prior          = sort(unique(o.uniqCond(:,1)),'descend');
cu             = o.uniqCond(:,2);      %coh
pu             = o.uniqCond(:,1);      %prior

%predicted densities by bin by 202 conditions
Ntrials        = 5000;
Ndots	       = 328;
Ns		       = 360;

%precompute stuffs for speed
pm0 = zeros(nCon,3);
pm0(cu==coh(1),2)   = simp(1);
pm0(cu==coh(2),2)   = simp(1);
pm0(cu==coh(3),2)   = simp(1);
pm0(pu==prior(1),1) = simp(2);
pm0(pu==prior(2),1) = simp(3);
pm0(pu==prior(3),1) = simp(4);
pm0(pu==prior(4),1) = simp(5);
kappa    = pm0(:,2);                  %sensory strength
kappa_s  = pm0(:,1);                  %prior strength
prand    = simp(6);                   %random estimate proba
km       = simp(7);                   %motor precision
s        = SLde2r(o.uniqCond(:,3),0); %s
su24     = unique(s(cu==coh(1)));     %dir. for coh 24
su12     = unique(s(cu==coh(2)));
su6      = unique(s(cu==coh(3)));
Nsu24    = length(su24);
Nsu12    = length(su12);
Nsu6     = length(su6);
x        = nan(Ntrials,Ndots,Ns);
coherent = nan(Ntrials,Ndots,Ns);
like     = nan(Ns,Ntrials);
svec     = SLde2r(1:1:360,0)';                               %hyp. directions
pmean    = SLde2r(225,0);                                           %high precision
motevec  = SLde2r(0:1:359,1);
[predEsttmp,pos]  = sort(SLde2r(1:1:360,1));      %estimate space
tailtmp   = (1 - coh)/2/pi;
amp1tmp   = coh/2/pi;
rotdeg24  = SLvectors2signedAngle(su24(1),su24,'radian');%rotation
rotdeg12  = SLvectors2signedAngle(su12(1),su12,'radian');
rotdeg6   = SLvectors2signedAngle(su6(1),su6,'radian');
rotrad24  = SLde2r(rotdeg24,1);
rotrad12  = SLde2r(rotdeg12,1);
rotrad6   = SLde2r(rotdeg6,1);
rotrad    = [{rotrad24} {rotrad12} {rotrad6}];
rotdeg    = [{rotdeg24} {rotdeg12} {rotdeg6}];
su        = [{su24} {su12} {su6}];
Nsu       = [{Nsu24} {Nsu12} {Nsu6}];
sPest     = nan(1,4*sum([Nsu{:}]));
prPest    = nan(1,4*sum([Nsu{:}]));
cohPest   = nan(1,4*sum([Nsu{:}]));
w{1}      = 1:4*Nsu{1};
w{2}      = 4*Nsu{1}+1 : 4*sum([Nsu{1:2}]);
w{3}      = 4*sum([Nsu{1:2}])+1 : 4*sum([Nsu{1:3}]);
M         = zeros(Ns,Ntrials);
svecStep  = SLvectors2signedAngle(svec(2),svec(1),'radian');
for i = 1 : length(coh)
    tmp = repmat([80 40 20 10],Nsu{i},1); %each prior (col)
    prPest(1,w{i}) = tmp(:);
    
    tmp = repmat(coh(i),1,4*Nsu{i});      %each coh (col)
    cohPest(1,w{i})= tmp;                 %coh
    
    tmp = repmat(su{i}',1,4);             %each direction (col)
    sPest(1,w{i}) = tmp;                  %stim. directions
    
    %calc. like density shifts
    shifttmp{i} = round(round(rotdeg{i})/svecStep);
    
    %warning
    if mod(round(rotdeg{i}),svecStep)
        fprintf('(slMakePredictionsWJM) Something"s wrong. Please change Ns. Ns=72 usually works fine.')
    end
    
end
[mcond,ix]= sortrows([prPest' cohPest' round(SLra2d(sPest))'],[1 2]);          %cond matr

%Match data & pred conditions
[~,nix]      = setdiff(mcond,o.uniqCond,'rows');  %non existing conditions out
Mcond        = mcond;
Mcond(nix,:) = [];
ix(nix)      = [];

%logl and predictions
%--------------------
tic
[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirlgTprior(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
o.duration = toc;
o.logl 	   = logl;

[~,bins]  = histc(0:1:360,xpdf);  %bins
bins(end) = [];
pPred     = nan(o.nCond,length(dataBins) - 1);
for i = 1 : o.nCond
    pPred(i,:) = SLcumSumInBin(PestGivModel(i,:),bins')';
end
pPred = pPred';


%backup
%-------
save('dataslMakePredictionsWJM')
fprintf('(dataslMakePredictionsWJM) Saving data in "dataslMakePredictionsWJM". \n')

%plot
%----
Bins = xpdf(2:end);
SLdrawModelsPredictionHistwithDataCentered(Pdata,Bins,pPred,db.stimFeatureDeg,sStrg,db.Pstd,priorModes,...
    o.uniqCond,'vonMisesPrior','experiment','vonMisesPrior');

%backup figures
myFig = SLgetFigures;
for i = 1: length(myFig)
    mynm = sprintf('%s%i','slMakePredictionsWJMlgTpriorFig0',i);
    saveas(myFig(i),[mynm,'.eps'],'epsc')
    saveas(myFig(i),[mynm,'.fig'],'fig')
end
fprintf('%s \n','-------------------------------- fig -----------------------------')
fprintf('%s \n','(slMakePredictionsWJMlgTpriorFig0) Your figures have been saved in')
fprintf(pwd);

!cp *eps ~/DropBox/
!cp *fig ~/DropBox/







