%slFitWJM.m
%
%
% author: steeve laquitaine
%   date: 150910
%purpose: fit WJM model to motion direction estimation data
%
%usage:
%
%       o = slFitWJM({'sub01'},[3 0.7 2.7 8 33 1e-2 32],'experiment','vonMisesPrior','cmaes','debug','dataPath','~/Dropbox/myDropbox/Project_withJustin/data/datapsychophy/proj01_priorstrength')
%       o = slFitWJM({'sub01'},[3 0.7 2.7 8 33 1e-2 32],'experiment','vonMisesPrior','cmaes','sherlock')
%
%
%Inputs:
%             'dataPath','~/': set data path
%                  'sherlock': cluster (default path is ~/data)
%
%options
%            'fminsearch': (or 'cmaes') fitting procedure
%                 'debug': minimal simulation for speed
%         
%
%Description:
%       - Motion direction likelihood is the integrated likelihood densities of the local dots directions
%       - Bayesian integration of Like and prior into posterior
%       - Estimates are posterior mean
%
%to copy fit results from Sherlock:
% >> scp -r steeve@sherlock.stanford.edu:data/fitWJMsub03 ~/Desktop/

%!scp -r steeve@sherlock.stanford.edu:data/cmaes/BayesNoCard/datafit10_vonMisesPrior_BayesNoCard_cmaes.mat ~/data/dataPsychophy/proj01_priorStrength/modelfit/AICcmaes/model_Bayes_MAP/


function o = slFitWJM(subjects,initp,varargin)

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
%input data path
if slIsInput(varargin,'dataPath')
    dataPath = varargin{find(strcmp(varargin,'dataPath'))+1};
    slPrintfStr('slFitWJM',['Gathering data in' dataPath])
elseif slIsInput(varargin,'sherlock')
    %datapath is in home directory on Sherlock cluster
    fprintf('%s \n','(slFitWJM) Gathering data ...')
    dataPath = '~/data';
end

cd(dataPath)
varg = [varargin,'dataPath',dataPath];    %experiment and path

%Organize data (estimates and task conditions)
db = SLMakedatabank(subjects,varg);     

%folder where fit results are backed up
mkdir(['fitWJM/fitWJM' subjects{1}])                      
cd(['fitWJM/fitWJM' subjects{1}])
fprintf('%s %s %s \n','(slFitWJM) Saving fit results in ', ['"data/fitWJM/fitWJM',subjects{1}],'"')

%collect data
data           = db.estimatesDeg;     %estimates
sStrg          = db.stimStrength;     %coherence
s              = db.stimFeatureDeg;   %direction
priorModes     = [db.priormodes{:}];  %prior modes
pstd           = db.Pstd;             %prior strength

%plot data densities by bin for each task condition (202)
dataBins       = 0:10:360;
[Pdata,xpdf,o] = slmakeDataDist(dataBins,db.estimatesDeg,db.stimFeatureDeg,db.stimStrength,db.Pstd,[db.priormodes{:}],'vonMisesPrior',varg);
o.Pdata        = Pdata;   %data density
o.posC         = o.posC;  %conditions
nCon           = o.nCond; %# of conditions
coh            = sort(unique(o.uniqCond(:,2)),'descend'); %coherence
prior          = sort(unique(o.uniqCond(:,1)),'descend'); %prior
cu             = o.uniqCond(:,2);      %unique coherence
pu             = o.uniqCond(:,1);      %unique prior

%predicted densities by bin for each task condition (202)
if slIsInput(varargin,'debug')
    Ntrials = 10;
    Ndots   = 10;
    Ns	   = 360;
else    
    Ntrials  = 5000;
    Ndots    = 328;
    Ns	     = 360;
end

%precompute stuffs for speed
pm0 = zeros(nCon,3);
pm0(cu==coh(1),2)   = initp(1);
pm0(cu==coh(2),2)   = initp(1);
pm0(cu==coh(3),2)   = initp(1);
pm0(pu==prior(1),1) = initp(2);
pm0(pu==prior(2),1) = initp(3);
pm0(pu==prior(3),1) = initp(4);
pm0(pu==prior(4),1) = initp(5);
kappa    = pm0(:,2);                  %sensory strength
kappa_s  = pm0(:,1);                  %prior strength
prand    = initp(6);                  %random estimate proba
km       = initp(7);                  %motor precision
s        = SLde2r(o.uniqCond(:,3),0); %s
su24     = unique(s(cu==coh(1)));     %unique directions for coh 24
su12     = unique(s(cu==coh(2)));     %unique directions for coh 12
su6      = unique(s(cu==coh(3)));     %unique directions for coh 6
Nsu24    = length(su24);
Nsu12    = length(su12);
Nsu6     = length(su6);
x        = nan(Ntrials,Ndots,Ns);
coherent = nan(Ntrials,Ndots,Ns);
like     = nan(Ns,Ntrials);
svec     = SLde2r(1:1:360,0)';                    %hyp. directions
pmean    = SLde2r(225,0);                         %prior mean
motevec  = SLde2r(0:1:359,1);
[predEsttmp,pos]  = sort(SLde2r(1:1:360,1));      %estimate space
tailtmp   = (1 - coh)/2/pi;
amp1tmp   = coh/2/pi;

%angle of rotation separating motion directions
rotdeg24  = SLvectors2signedAngle(su24(1),su24,'radian');
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
        fprintf('(slFitWJM) Something"s wrong. Please change Ns. Ns=72 usually works fine.')
    end
    
end
[mcond,ix]= sortrows([prPest' cohPest' round(SLra2d(sPest))'],[1 2]);%cond matr

%Match data & pred conditions
[~,nix]      = setdiff(mcond,o.uniqCond,'rows'); %non existing conditions out
Mcond        = mcond;
Mcond(nix,:) = [];
ix(nix)      = [];

%---
%fit
%---
%init fit backup variables
history.params    = [];
history.iteration = [];
history.funccount = [];
history.fval      = [];
tic

%fminsearch
if slIsInput(varargin,'fminsearch')
    slPrintfStr('slFitWJM','Fitting data with fminsearch procedure')
    options = optimset('OutputFcn', @myoutput,'MaxIter',2000,'MaxFunEvals',10*2000,'Display','iter');
    %fit
    [fitP,neglogl,exitflag,output] = fminsearch(@(fitP) ...
                                                SLgetLogL_WJMrotatedDir(data,s,cu,kappa,kappa_s,prand,km,fitP,o.posC,o.uniqCond),initp,options);
end
%cmaes
if slIsInput(varargin,'cmaes')
    slPrintfStr('slFitWJM','Fitting data with cmaes procedure')
    
    %options
    options         = cmaes;
    options.LBounds = 0; %no < 0
    options.TolFun  = 1e-4; %stop when obj(t2) - obj(t1) < 1e-4
    
    %debug mode
    if slIsInput(varargin,'debug')
        options.MaxIter = 1;
    end
    
    %have to set those off otherwise parfor does not work
    options.DispFinal     = 0;
    options.DispModulo    = 0;
    options.SaveVariables = 0;
    options.SaveFilename  = 0;
    options.LogModulo     = 0;
    options.LogTime       = 0;
    options.LogFilenamePrefix = 0;
    options.LogPlot       = 0;

    %cmaes
    %    [fitP,neglogl,exitflag,output] = cmaes(@(fitP) ...
    %                                           SLgetLogL_WJMrotatedDirCmaes(fitP,data,s,c,kappa,kappa_s,prand,km,posC,Cond))
    [fitP,neglogl,counteval,exitflag,output] = cmaes('SLgetLogL_WJMrotatedDirCmaes',initp,sqrt(var(initp)),options,data,s,cu,kappa,kappa_s,prand,km,o.posC,o.uniqCond);

end


o.duration  = toc;
o.fitP      = fitP;
o.neglogl   = neglogl;
o.logl 	    = - neglogl;
o.exitflag  = exitflag;
o.outputFit = output;

%best predictions to plot
[~,~,PestGivModel] = SLgetLogL_WJMrotatedDir(data,s,cu,kappa,kappa_s,prand,km,o.fitP,o.posC,o.uniqCond);
o.PpredRaw = PestGivModel;
[~,bins]   = histc(0:1:360,xpdf);  %bins
bins(end)  = [];
pPred      = nan(o.nCond,length(dataBins) - 1);
for i = 1 : o.nCond
    pPred(i,:) = SLcumSumInBin(PestGivModel(i,:),bins')';
end
pPred = pPred';
o.PpredBinned = pPred;

%backup
%------
exp = varargin{find(sum(strcmp(varargin,'experiment')))+1};
save(['datafit', subjects{1}(end-1:end),'_',exp,'_WJM'],'o')
fprintf('%s %s \n','(slFitWJM) Saving data as "',['datafit' subjects{1}(end-1:end) '_' exp '_WJM "'])

%plot
%----
%Bins = xpdf(2:end);
%o.Bins = Bins;
%SLdrawModelsPredictionHistwithDataCentered(Pdata,Bins,pPred,db.stimFeatureDeg,sStrg,db.Pstd,priorModes,...
%    o.uniqCond,'vonMisesPrior','experiment','vonMisesPrior');

    %fit output
    function stop = myoutput(x,optimvalues,state)

        stop = false;
        if isequal(state,'iter')
            history.params    = [history.params    ; x                    ];
            history.iteration = [history.iteration ; optimvalues.iteration];
            history.funccount = [history.funccount ; optimvalues.funccount];
            history.fval      = [history.fval      ; optimvalues.fval     ];
        end
        
        %save filename fname
        save('fitHistory','history')
    end
end
