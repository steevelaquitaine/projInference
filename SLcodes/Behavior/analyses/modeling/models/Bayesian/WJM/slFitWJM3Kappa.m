%slFitWJM3Kappa.m
%
%
% author: steeve laquitaine
%   date: 150910
%purpose: fit WJM model to motion direction estimation data
%
%usage:
%
%       o = slFitWJM3Kappa({'sub01'},[40 5 5 2.77 4.5 15 300 1e-2 16],'experiment','vonMisesPrior')

function o = slFitWJM3Kappa(subjects,initp,varargin)

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
% mfilename    = SLgetActivemFile;
% myMfile      = SLbackupMfileInWSpace(mfilename);
% o.mfilename  = mfilename;
% o.myMfile    = myMfile;

%--------
%get data
%--------
fprintf('%s \n','(slFitWJM) Gathering data ...')
%fprintf('%s \n','(slMakePredictionsWJM) Please set dataPath ...')
dataPath     = '~/';
cd(dataPath)
varg 		 = [varargin,'dataPath',dataPath];    %experiment and path
db           = SLMakedatabank(subjects,varg);
mkdir(['fitWJM' subjects{1}])
cd(['fitWJM' subjects{1}])
fprintf('%s %s %s \n','(slFitWJM) Saving fit results in', ['"data/fitWJM',subjects{1}],'"')
data           = db.estimatesDeg;     %estimates
sStrg          = db.stimStrength;     %coherence
s              = db.stimFeatureDeg;   %direction
priorModes     = [db.priormodes{:}];  %prior modes
pstd           = db.Pstd;             %prior strength

%data densities by bin by 202 conditions
dataBins       = 0:10:360;
[Pdata,xpdf,o] = slmakeDataDist(dataBins,db.estimatesDeg,db.stimFeatureDeg,db.stimStrength,db.Pstd,[db.priormodes{:}],'vonMisesPrior');
o.Pdata        = Pdata;
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
pm0(cu==coh(1),2)   = initp(1);
pm0(cu==coh(2),2)   = initp(2);
pm0(cu==coh(3),2)   = initp(3);
pm0(pu==prior(1),1) = initp(4);
pm0(pu==prior(2),1) = initp(5);
pm0(pu==prior(3),1) = initp(6);
pm0(pu==prior(4),1) = initp(7);
kappa    = pm0(:,2);                  %sensory strength
kappa_s  = pm0(:,1);                  %prior strength
prand    = initp(8);                   %random estimate proba
km       = initp(9);                   %motor precision
initp     = initp;
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
svec     = (0:2*pi/Ns:2*pi-2*pi/Ns)';       %hyp. directions
pmean    = 3.93;
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
        fprintf('(slFitWJM) Something"s wrong. Please change Ns. Ns=72 usually works fine.')
    end
    
end
[mcond,ix]= sortrows([prPest' cohPest' round(SLra2d(sPest))'],[1 2]);          %cond matr

%Match data & pred conditions
[~,nix]      = setdiff(mcond,o.uniqCond,'rows');  %non existing conditions out
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

%options = optimset('OutputFcn', @myoutput,'MaxIter',1,'MaxFunEvals',10,'Display','iter'); %quick fit
options = optimset('OutputFcn', @myoutput,'MaxIter',2000,'MaxFunEvals',10*2000,'Display','iter');

tic
[fitP,neglogl,exitflag,output] = fminsearch(@(fitP) SLgetLogL_WJMrotatedDir3Kappa(data,s,cu,kappa,kappa_s,prand,km,fitP,o.posC,o.uniqCond),initp,options);
o.duration  = toc;
o.fitP      = fitP;
o.neglogl 	= neglogl;
o.logl 	    = - neglogl;
o.exitflag  = exitflag;
o.outputFit = output;

%best pred
[~,~,PestGivModel] = SLgetLogL_WJMrotatedDir(data,s,cu,kappa,kappa_s,prand,km,fitP,o.posC,o.uniqCond);
o.PpredRaw= PestGivModel;
[~,bins]  = histc(0:1:360,xpdf);  %bins
bins(end) = [];
pPred     = nan(o.nCond,length(dataBins) - 1);
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
Bins = xpdf(2:end);
o.Bins = Bins;
SLdrawModelsPredictionHistwithDataCentered(Pdata,Bins,pPred,db.stimFeatureDeg,sStrg,db.Pstd,priorModes,...
    o.uniqCond,'vonMisesPrior','experiment','vonMisesPrior');

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
