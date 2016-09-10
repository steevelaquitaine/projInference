
%slMakePredictionsWJM0.m
%
%
% author: steeve laquitaine
%   date: 150724
%purpose: calculate motion direction estimates predictions with WJma
%         model and the logl of data given the model.
%
%  usage: 
%
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJM0({'sub01'},[60 20 1 .8 2.77 8.7 33.3 1e-6 1000],'experiment','vonMisesPrior')
%	   [Pdata,pPred,Bins,o] = slMakePredictionsWJM0({'sub01'},[40 5 5 2.77 4.5 15 300 1e-2 16],'experiment','vonMisesPrior')
%	   [Pdata,pPred,Bins,o] = slMakePredictionsWJM0({'sub01'},[40 40 40 0 0 0 0 1e-2 16],'experiment','vonMisesPrior')
%	   
%      
%
%10 cores:
%duration: 5 min (5000 trials)
%duration: 3 min (2500 trials)
%duration: 2 min (1250 trials)

function [Pdata,pPred,Bins,o] = slMakePredictionsWJM0(subjects,simp,varargin)

tic
% matlabpool 11	

%--------------
%backup .m file 
%--------------
%mfilename = SLgetActivemFile;
%myMfile = SLbackupMfileInWSpace(mfilename);
%o.mfilename  = mfilename;
%o.myMfile = myMfile;

%-----------------
%set path and data
%-----------------
fprintf('%s \n','(slMakePredictionsWJM) Gathering data ...')
fprintf('%s \n','(slMakePredictionsWJM) Please set dataPath ...')
%dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');
%dataPath = '~/data/dataPsychophy/proj01_priorStrength';

%%%%%%%%% to uncomment %%%%%%%%% 

%dataPath = '~/Dropbox/myDropbox/Project_withJustin/data/dataPsychophy/proj01_priorStrength';
%cd(dataPath)
%varg 		 = [varargin,'dataPath',dataPath];    %experiment and path
%db           = SLMakedatabank(subjects,varg);

%%%%%%%%% to uncomment %%%%%%%%% 

%%%%%%%%% load data fast (for debugging) %%%%%%%%% 

load datbank
db = databank;

%%%%%%%%% load data fast (for debugging) %%%%%%%%% 


data         = db.estimatesDeg;           		  %subjects estimates
%pstd         = db.Pstd;                   		  %task priors std
sStrg        = db.stimStrength;           		  %stimulus coherence
s            = db.stimFeatureDeg;                 %stimulus direction
priorModes   = [db.priormodes{:}];        		  %prior modes


%-----------------
%data distribution
%-----------------
dataBins  = 0 : 10 : 360;		  %bins of data densities

%data densities 202 conditions
[Pdata,xpdf,o] = slmakeDataDist(dataBins,data,s,sStrg,db.Pstd,priorModes,'vonMisesPrior');            

%p(data) by bin by 202 conditions
pm0 = o.uniqCond;
pm0(o.uniqCond==0.24) = simp(1);
pm0(o.uniqCond==0.12) = simp(2);
pm0(o.uniqCond==0.06) = simp(3);
pm0(o.uniqCond==80)   = simp(4);
pm0(o.uniqCond==40)   = simp(5);
pm0(o.uniqCond==20)   = simp(6);
pm0(o.uniqCond==10)   = simp(7);
kappa   = pm0(:,2);                  %sensory strength
kappa_s = pm0(:,1);                  %prior strength
prand   = simp(8);                   %random estimate proba
km      = simp(9);                   %motor precision
s       = SLde2r(o.uniqCond(:,3),0); %s
c       = o.uniqCond(:,2);           %coh

%--------------
%Preallocation
%--------------
global predEst     %estimate space
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
Ntrials = 250;
Ndots	= 328;
Ns		= 100;
su24  = unique(s(c==0.24));    %dir. for coh 24
su12  = unique(s(c==0.12));    %dir. for coh 24
su6   = unique(s(c==0.06));    %dir. for coh 24
Nsu24 = length(su24);
Nsu12 = length(su12);
Nsu6  = length(su6);
predEst  = 1 : 1 : 360;
x        = nan(Ntrials,Ndots,Ns);    
coherent = nan(Ntrials,Ndots,Ns); 
like     = nan(Ns,Ntrials);  


%-----------
%predictions
%-----------
%[logl,~,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,simp,o.posC);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDir(data,s,c,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest0(data,s,c,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
tocs = toc;
o.duration = tocs;

o.logl 	  = logl;
[~,bins]  = histc(0 : 1 : 360, xpdf);  %make bins
bins(end) = [];
bins      = bins'; %space for pred and data densities
pPred     = nan(o.nCond,length(dataBins) - 1);

for i = 1 : o.nCond
    pPred(i,:) = SLcumSumInBin(PestGivModel(i,:),bins)';
end
pPred = pPred';

save('dataslMakePredictionsWJM')
fprintf('(dataslMakePredictionsWJM) Saving data in "dataslMakePredictionsWJM". \n')


%----
%plot
%----
Bins = xpdf(2:end);
SLdrawModelsPredictionHistwithDataCentered0(Pdata,Bins,pPred,db.stimFeatureDeg,sStrg,db.Pstd,priorModes,...
   o.uniqCond,'vonMisesPrior','experiment','vonMisesPrior')