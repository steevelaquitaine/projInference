
%slMakePredictionsWJM5.m
%
%
% author: steeve laquitaine
%   date: 150724
%purpose: calculate motion direction estimates predictions with WJma
%         model and the logl of data given the model.
%
%  usage: 
%
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJM5({'sub01'},[60 20 1 .8 2.77 8.7 33.3 1e-6 1000],'experiment','vonMisesPrior')
%	   [Pdata,pPred,Bins,o] = slMakePredictionsWJM5({'sub01'},[40 5 5 2.77 4.5 15 300 1e-2 16],'experiment','vonMisesPrior')
%	   [Pdata,pPred,Bins,o] = slMakePredictionsWJM5({'sub01'},[40 40 40 0 0 0 0 1e-2 16],'experiment','vonMisesPrior')
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJM5({'sub01'},[300 300 300 0 0 0 0 1e-2 16],'experiment','vonMisesPrior')
%	   
%      
%
%10 cores:
%duration: 5 min (5000 trials)
%duration: 3 min (2500 trials)
%duration: 2 min (1250 trials)

function [Pdata,pPred,Bins,o] = slMakePredictionsWJM5(subjects,simp,varargin)

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
global svect3D
tic

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
% dataPath = '~/Dropbox/myDropbox/Project_withJustin/data/dataPsychophy/proj01_priorStrength';
% cd(dataPath)
% varg 		 = [varargin,'dataPath',dataPath];    %experiment and path
% db           = SLMakedatabank(subjects,varg);
%%%%%%%%% load data fast (for debugging) %%%%%%%%% 
load datbank
db = databank;
%%%%%%%%% load data fast (for debugging) %%%%%%%%% 

data                = single(db.estimatesDeg);     %subjects estimates
sStrg               = single(db.stimStrength);     %stimulus coherence
s                   = single(db.stimFeatureDeg);   %stimulus direction
priorModes          = single([db.priormodes{:}]);  %prior modes
pstd                = single(db.Pstd);             %prior strength

%data densities by bin by 202 conditions
dataBins            = single(0:10:360);  %data densities bin
[Pdata,xpdf,o]      = slmakeDataDist(dataBins,db.estimatesDeg,db.stimFeatureDeg,db.stimStrength,db.Pstd,[db.priormodes{:}],'vonMisesPrior');            
o.posC              = single(o.posC);
nCon                = single(o.nCond);
coh                 = single(sort(unique(o.uniqCond(:,2)),'descend'));
prior               = single(sort(unique(o.uniqCond(:,1)),'descend'));
cu                  = single(o.uniqCond(:,2));          %coh
pu                  = single(o.uniqCond(:,1));           %prior

%predicted densities by bin by 202 conditions
Ntrials             = single(5000);
Ndots               = single(328);
Ns                  = 100;
pm0                 = zeros(nCon,3,'single');
pm0(cu==coh(1),2)   = simp(1);
pm0(cu==coh(2),2)   = simp(2);
pm0(cu==coh(3),2)   = simp(3);
pm0(pu==prior(1),1) = simp(4);
pm0(pu==prior(2),1) = simp(5);
pm0(pu==prior(3),1) = simp(6);
pm0(pu==prior(4),1) = simp(7);
kappa               = pm0(:,2);                  %sensory strength
kappa_s             = pm0(:,1);                  %prior strength
prand               = simp(8);                   %random estimate proba
km                  = simp(9);                   %motor precision
s                   = SLde2r(o.uniqCond(:,3),0); %s
su24                = unique(s(cu==coh(1)));     %dir. for coh 24
su12                = unique(s(cu==coh(2)));
su6                 = unique(s(cu==coh(3)));
Nsu24               = length(su24);
Nsu12               = length(su12);
Nsu6                = length(su6);
x                   = nan(Ntrials,Ndots,Ns);    
coherent            = nan(Ntrials,Ndots,Ns); 
like                = nan(Ns,Ntrials);  
svec                = single(0:2*pi/Ns:2*pi-2*pi/Ns)';        %hyp. directions
pmean               = single(3.93);                           %prior mean
svect3D             = double(repmat(svec,[1 Ntrials Ndots])); %Ns by Ntrials by Ndots

%predictions
%-----------
%[logl,~,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,simp,o.posC);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDir(data,s,c,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest2(data,s,c,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest3(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest5(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);


[~,bins]  = histc(0:1:360,xpdf);  %bins
bins(end) = [];
pPred     = nan(o.nCond,length(dataBins) - 1);
for i = 1 : o.nCond
    pPred(i,:) = SLcumSumInBin(PestGivModel(i,:),bins')';
end
pPred = pPred';


%backup
%-------
o.duration = toc;
o.logl 	  = logl;
save('dataslMakePredictionsWJM')
fprintf('(dataslMakePredictionsWJM) Saving data in "dataslMakePredictionsWJM". \n')


%plot
%----
Bins = xpdf(2:end);
SLdrawModelsPredictionHistwithDataCentered(single(Pdata),single(Bins),single(pPred),single(db.stimFeatureDeg),single(sStrg),single(db.Pstd),single(priorModes),...
 single(o.uniqCond),'vonMisesPrior','experiment','vonMisesPrior');

