

%slfitMaxLLH.m
%
%
% author: steeve laquitaine
%purpose: fit Bayesian model with maximum likelihood
%
%  usage:
%
%      slfitMaxLLH(database,fitP,initp,'vonMisesPrior',filename,output,varargin);
%
%Description:
%
%   algorithm : fminsearch, Nelder Mead
%
%nested functions : slfitLearningBayes, slgetloglLearningBayesianModel

function [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
    slfitMaxLLH(database,fitPinput,initp,filename,output,varargin)

fprintf('%s \n','(slfitMaxLLH) Fitting trial-data with max likelihood \n')

%case test data
if sum(strcmp(varargin{1},'testFakedata'))
    [data, disp, coh, pstd, priorModes] = setTestData;
end

%extract from database
fprintf('%s \n','(slfitMaxLLH) Extracting data and variables \n')
[output,data,disp,StimStrength,pstd,priorModes,trials] = extractvarFromDatabase(database);

%fit
if isempty(fitPinput)==1
    
    %set initial parameters
    %InitParams_freePriors
    %(1) Input parameters (best matching parameters)
    %Input initial fit parameters (typically best matching initial values)
    %check that all initial parameters are input
    if length(initp)~=12
        fprintf('%s \n',['(slfitMaxLLH) Initial parameters are missing...',...
            'parameters for learning Bayesian model is:', '\n 3 kllh', '\n 4 klearnt','\n 1 kcardinal',...
            '\n 1 fractRand','\n 1 motorNoise','\n 1 weighTail','\n 1 learning tau'])
        keyboard
    end
    kllh       = initp(1:3);
    klearnt    = initp(4:7);
    kcardinal  = initp(8);
    fractRand  = initp(9);
    motorNoise = initp(10);
    weightTail = initp(11);
    tau        = initp(12);
    k0         = [kllh klearnt kcardinal fractRand motorNoise weightTail tau];
    output.initp = k0;

    %no cardinal prior
    if isnan(kcardinal); TheModel = 'withoutCardinal'; else TheModel = 'withCardinal'; end    
    
    fitPbkp = nan(size(k0,1),size(k0,2));
    negLoglbkp = nan(size(k0,1),1);
    fprintf('%s \n','(slfitMaxLLH) Algo : fminsearch \n')
    options = optimset('MaxIter',200*size(k0,1),'MaxFunEvals',10*200*size(k0,1));    
    
    %print expected duration
    duration = 19;
    %expDuration =  duration*10*200*9/3600;
    expDuration =  duration*100/3600;
    fprintf('%s %.2f %s /n','Expected duration: ', expDuration,'hours')
    
    [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp,fithistory] = slfitLearningBayes(output.initp,data,disp,StimStrength,...
        pstd,priorModes,trials,TheModel,options,varargin{:});
    output.freePriors = 1;
    output.fixedPriors = 0;    
    
    %backup
    output.fitPbkp = fitPbkp;
    output.negLoglbkp = negLoglbkp;
    output.outputFitbkp = outputFitbkp;
    output.exitflagbkp  = exitflagbkp;
    output.fithistory = fithistory;
    
    %Max likelihood. To have an idea of what log likelihood we should
    %obtain our reasoning is that if one cannot predict at all motion
    %direction at a given trial, the probability that any given direction
    %be chosen is 1/360. (we consider estimate with 1 degree resolution).
    %So loglikelihood should be log(1/360)=-5.88. This mutiplied (because
    %of log) by our ~5500 trials, we get negLogl~32000. Any ability to
    %somewhat predict better motion direction better than chance yields
    %soemthing lower than 32000.
    [negLogl,position] = min(negLoglbkp(:));
    i = ind2sub(size(negLoglbkp),position);
    
    %optimal parameters and Logl_pertrial
    fitP.p = fitPbkp(i,:);
    if any(strcmp(varargin{:},'fminsearch'))
        [~,~,Logl_pertrialBestfit] = slgetloglLearningBayesianModel(data,disp,StimStrength,pstd,fitP.p,...
        priorModes,trials,TheModel,varargin{1});
    end
    
    [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp,fithistory] = slfitLearningBayes(initp,data,disp,StimStrength,...
        pstd,priorModes,trials,TheModel,options,varargin);
    
    %backup
    output.Logl_pertrialBestfit = Logl_pertrialBestfit;
    output.fitP = fitP.p;
    
    fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','klearnt80','klearnt40','klearnt20',...
        'klearnt10','Kcardinal','Prand','km'};
    
    %no cardinal prior
    if isnan(kcardinal); fitP.p(strcmp(fitP.nm,'kcardinal')) = NaN; end
    
    %simple search
    output.TrialPred = [];
    output.stdPa  = [];
    
    %save
    save(filename)
    
    %get logl for input fitp
elseif isempty(fitPinput)==0
    
    %no cardinal prior
    if isnan(fitPinput(8)); TheModel = 'withoutCardinal'; else TheModel = 'withCardinal'; end
    [negLogl,~,Logl_pertrialBestfit] = slgetloglLearningBayesianModel(data,disp,StimStrength,pstd,fitPinput,...
        priorShape,priorModes,TheModel,varargin{1});
    negLoglbkp = negLogl;
    fitPbkp = fitPinput;
    fitP.p = fitPinput;
    fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','klearnt80','klearnt40','klearnt20',...
        'klearnt10','Kcardinal','Prand','km'};
    %save
    save(filename)
end

%get test data
function [data, disp, coh, pstd, priorModes] = setTestData

fprintf('%s \n','---------------------------------------------------- \n')
fprintf('%s \n','(slfitMaxLLH) loading simulated data to check fitting \n')
fprintf('%s \n','---------------------------------------------------- \n')
data = evalin('base','pred');
disp = evalin('base','d');
coh  = evalin('base','coh');
pstd = evalin('base','pstd');
priorModes = evalin('base','priorModes');
%extract variables
function [output,data,disp,StimStrength,pstd,priorModes,trials] = extractvarFromDatabase(database)

%subjects' data and variables
data         = round(cell2mat(database.data(:,(strcmp(database.nm,'estimatedFeature'))==1)));
disp         = cell2mat(database.data(:,(strcmp(database.nm,'FeatureSample'))==1));
StimStrength = cell2mat(database.data(:,(strcmp(database.nm,'StimStrength'))==1));
pstd         = cell2mat(database.data(:,(strcmp(database.nm,'Pstd'))==1));
priorModes   = cell2mat(database.data(:,(strcmp(database.nm,'priormodes'))==1));
trials       = cell2mat(database.data(:,(strcmp(database.nm,'Trials'))==1));

%make 360 and 0 degrees same
data(data==0) = 360;

%drop missing data
data          = data(isnan(data)==0);
disp          = disp(isnan(data)==0);
StimStrength  = StimStrength(isnan(data)==0);
pstd          = pstd(isnan(data)==0);
priorModes    = priorModes(isnan(data)==0,:);
trials        = trials(isnan(data)==0,:);

%store
output.data = data;
output.disp = disp;
output.StimStrength = StimStrength;
output.pstd = pstd;
output.priorModes = priorModes;
output.trials = trials;
