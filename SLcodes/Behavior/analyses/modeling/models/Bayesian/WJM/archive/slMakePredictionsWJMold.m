
%slMakePredictionsWJMold.m
%
%
% author: steeve laquitaine
%   date: 150724
%purpose: calculate motion direction estimates predictions with WJma
%         model and the logl of data given the model.
%
%  usage: 
%
%      [Pdata,pPred,Bins,o] = slMakePredictionsWJMold({'sub01'},[60 20 1 .8 2.77 8.7 33.3 1e-6 1000],'experiment','vonMisesPrior')
%	   [Pdata,pPred,Bins,o] = slMakePredictionsWJMold({'sub01'},[40 5 5 2.77 4.5 15 300 1e-2 16],'experiment','vonMisesPrior')
%	  
%      
%
%10 cores:
%duration: 5 min (5000 trials)
%duration: 3 min (2500 trials)
%duration: 2 min (1250 trials)

function [Pdata,pPred,Bins,o] = slMakePredictionsWJMold(subjects,simp,varargin)

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

%--------------
%backup .m file 
%--------------
%mfilename = SLgetActivemFile;
%myMfile = SLbackupMfileInWSpace(mfilename);
%o.mfilename  = mfilename;
%o.myMfile = myMfile;

%--------
%get data
%--------
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

data           = single(db.estimatesDeg);     %subjects estimates
sStrg          = single(db.stimStrength);     %stimulus coherence
s              = single(db.stimFeatureDeg);   %stimulus direction
priorModes     = single([db.priormodes{:}]);  %prior modes
pstd           = single(db.Pstd);             %prior strength

%data densities by bin by 202 conditions
dataBins       = single(0:10:360);  %data densities bin
[Pdata,xpdf,o] = slmakeDataDist(dataBins,db.estimatesDeg,db.stimFeatureDeg,db.stimStrength,db.Pstd,[db.priormodes{:}],'vonMisesPrior');            
o.posC         = single(o.posC);
nCon           = single(o.nCond);
coh            = single(sort(unique(o.uniqCond(:,2)),'descend'));
prior          = single(sort(unique(o.uniqCond(:,1)),'descend'));
cu             = single(o.uniqCond(:,2));          %coh
pu             = single(o.uniqCond(:,1));           %prior

%predicted densities by bin by 202 conditions
Ntrials  = single(5000);
Ndots	 = single(328);
Ns		 = 360;%100; for SLgetLogL_WJMrotatedDirTest7 only
%Ns		 = 100;      for others (higher lead to memory load)

%precompute stuffs for speed
pm0 = zeros(nCon,3,'single');
pm0(cu==coh(1),2)   = simp(1);
pm0(cu==coh(2),2)   = simp(2);
pm0(cu==coh(3),2)   = simp(3);
pm0(pu==prior(1),1) = simp(4);
pm0(pu==prior(2),1) = simp(5);
pm0(pu==prior(3),1) = simp(6);
pm0(pu==prior(4),1) = simp(7);
kappa    = pm0(:,2);                  %sensory strength
kappa_s  = pm0(:,1);                  %prior strength
prand    = simp(8);                   %random estimate proba
km       = simp(9);                   %motor precision
simp     = double(simp);
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
svec     = double(0:2*pi/Ns:2*pi-2*pi/Ns)';       %hyp. directions
pmean    = single(3.93);
motevec  = SLde2r(0:1:359,1);
[predEsttmp,pos]  = sort(SLde2r(1:1:360,1));      %estimate space
tailtmp   = (1 - double(coh))/2/pi;
amp1tmp   = double(coh)/2/pi;
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
M         = zeros(Ns,Ntrials,'single');
svecStep  = SLvectors2signedAngle(svec(2),svec(1),'radian'); 
for i = 1 : length(coh)
    tmp = repmat([80 40 20 10],Nsu{i},1); %each prior (col)
    prPest(1,w{i}) = tmp(:);
    
    tmp = repmat(coh(i),1,4*Nsu{i}); %each coh (col)
    cohPest(1,w{i})= tmp;            %coh
    
    tmp = repmat(su{i}',1,4);        %each direction (col)
    sPest(1,w{i}) = tmp;             %stim. directions
    
    %calc. like density shifts    
    shifttmp{i} = round(round(rotdeg{i})/svecStep);
    
    %warning
    if mod(round(rotdeg{i}),svecStep)
        fprintf('(slMakePredictionsWJM) Something"s wrong. Please change Ns. Ns=72 usually works fine.')
    end

end
[mcond,ix]= sortrows([prPest' cohPest' round(SLra2d(sPest))'],[1 2]);          %cond matr

%Match data and model conditions
[~,nix]   = setdiff(single(mcond),single(o.uniqCond),'rows');  %non existing conditions out
Mcond     = mcond;
Mcond(nix,:) = [];                 
ix(nix)      = [];

%predictions
%-----------
tic
%[logl,~,PestGivModel] = SLgetLogL_WJM(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDir(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest2(data,s,c,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest3(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest4(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
%[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest6(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
[logl,~,PestGivModel] = SLgetLogL_WJMrotatedDirTest7(data,s,cu,kappa,kappa_s,prand,km,simp,o.posC,o.uniqCond);
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
SLdrawModelsPredictionHistwithDataCentered(single(Pdata),single(Bins),single(pPred),single(db.stimFeatureDeg),single(sStrg),single(db.Pstd),single(priorModes),...
 single(o.uniqCond),'vonMisesPrior','experiment','vonMisesPrior');

