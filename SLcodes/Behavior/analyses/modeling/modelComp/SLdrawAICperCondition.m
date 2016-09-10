
% SLdrawAICperCondition.m
%
%
%
%     author: steeve laquitaine
%       date: 140526 (last update 150105)
%    purpose: draw AIC or cross-validated R2 for comparison of different
%             model performances for different subjects, for different
%             experiments, for different experimental conditions.
%
%description:
%
%             AIC = 2 * (numfitPthisModel - loglall)
%
%   with loglall = logl(datat1/model1) + logl(datat2/model1) + ....
%
%
%usage:
%
%Selected model AIC based on entire dataset.
%
% SLdrawAICperCondition({'sub01',...
%           'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},...
%            'experiment','vonMisesPrior','models',...
%            {'model_CompDiv',... 
%            'model_Bayes_MAP',... 
%            'model_Bayes_Sampling',...
%            'model_Bayes_MAP_withCard',...
%            'model_Bayes_Sampling_withCard',... 
%            'model_Bayes_WJM',...
%            'model_Bayes_MAP_FatTailPrior',...
%            'model_Bayes_WJMtailedPrior'},'AIC','DeltaI','subjectAverage');
%
%
%cmaes fit so far
% SLdrawAICperCondition({'sub01',...
%           'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'},...
%            'experiment','vonMisesPrior','models',...
%            {'model_CompDiv',... 
%            'model_Bayes_MAP',... 
%            'model_Bayes_WJM',...
%            'model_Bayes_Sampling'},'AICcmaes','DeltaI','subjectAverage');
% 
%
%list of analyses
%----------------
%
%         AIC: draw the AIC obtained from maximum likelihood fitting of
%                               different models with the data.
%
%         CrossValidated R^2: draw the mean cross validated R2 for each
%                               subject,model obtain from least square fit 
%                               of the data
%
%varargs options:
%----------------
%     'dataPath' followed by path
%
%     'experiment': 'vonMisesPrior': circular prior experiment
%     'experiment': 'bimodalPrior': mixture of circular prior experiment
%     'models': followed by 
%                     'model_Bayes_MAP'
%                     'model_Bayes_Sampling'
%                     'model_Bayes_Sampling_withCard'
%                     'model_Bayes_MAP_FatTailPrior'
%                     'model_Bayes_Sampling_FatTailPriorAndLLH'
%                     'model_Bayes_Sampling_FatTailPrior_ChWT_NoCard'
%                     'model_Bayes_MAP_FatTailPrior_withCard'
%                     'model_Bayes_MAP_withCard'
%                     'model_CompDiv'
%                     'model_CompDiv_withCard'
%                     'model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails'
%                     'model_Bayes_WJM'
%                     'model_Bayes_WJMtailedPrior'
%
%     'loglEachCondition': plot the log likelihood of the data given the
%                          model for each experimental conditions.
%                   'AIC': model comparison with AIC. Lowest AIC model fits data best.
%             'logDeltaI': plot log DeltaI
%                'DeltaI': plot DeltaI
%            'CrossValR2': model comparison with Cross-validated R^2. Highest R^2
%                          explains data best.
%       'CrossValmaxlogl': model comparison with Cross-validated maxlogl.
%                          Highest maxlogl explains data best.
%    'distanceToPrior=x': analyse only trials where displayed direction is
%                          x deg (e.g., 90) from prior mean
%
%
%
%USAGE:
%
%
%
%AIC (log10 scaled) -  subject average
%%-------------------------------------
% pathLab =  '/Users/steeve/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2';
% pathHome = '/Users/steeve_laquitaine/Dropbox/Project_withJustin/data/dataPsychophy/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior';
% pathHome = '/Users/steeve_laquitaine/Dropbox/Project_withJustin/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/';
%
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01'},...
%     'experiment','vonMisesPrior',...
%     'models',{'model_CompDiv','model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails'},...
%     'AIC','distanceToPrior=90');
% 
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub06'},...
%     'experiment','vonMisesPrior',...
%     'models',{'model_CompDiv',...
%     'model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails'},...
%     'AIC','distanceToPrior=90');
% 
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10'...
%     'sub11','sub12'},...
%     'experiment','vonMisesPrior',...
%     'dataPath',pathLab,...
%     'subjectAverage',...
%     'models',{'model_CompDiv',...
%     'model_Bayes_WJM',...
%     'model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails',...
%     'model_CompDiv_withCard',...
%     'model_Bayes_MAP',...
%     'model_Bayes_Sampling',...
%     'model_Bayes_Sampling_withCard',...
%     'model_Bayes_MAP_FatTailPrior',...
%     'model_Bayes_MAP_withCard'},'AIC','distanceToPrior=90');
%
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10'...
%     'sub11','sub12'},...
%     'experiment','vonMisesPrior',...
 %    'models',{'model_CompDiv',...
 %    'model_Bayes_WJM',...
 %    'model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails',...
 %    'model_CompDiv_withCard',...
 %    'model_Bayes_MAP',...
 %    'model_Bayes_Sampling',...
 %    'model_Bayes_Sampling_withCard',...
 %    'model_Bayes_MAP_FatTailPrior',...
 %    'model_Bayes_Sampling_FatTailPriorAndLLH',...
 %    'model_Bayes_Sampling_FatTailPrior_ChWT_NoCard',...
 %    'model_Bayes_MAP_FatTailPrior_withCard',...
 %    'model_Bayes_MAP_withCard'},'AIC','distanceToPrior=0');
%
 %[loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
 %    'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10'...
 %    'sub11','sub12'},...
 %    'experiment','vonMisesPrior',...
 %     'subjectAverage',...
 %   'models',{'model_CompDiv',...
 %   'model_Bayes_WJM',...
 %   'model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails',...
 %   'model_CompDiv_withCard',...
 %   'model_Bayes_MAP',...
 %   'model_Bayes_Sampling',...
 %   'model_Bayes_Sampling_withCard',...
 %   'model_Bayes_MAP_FatTailPrior',...
 %   'model_Bayes_MAP_withCard'},'AIC','distanceToPrior=0');
%
% 
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub02'},...
%     'experiment','vonMisesPrior',...
%     'models',{'model_Bayes_Sampling_withCard'},'AIC','distanceToPrior=0');
% 
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10'...
%     'sub11','sub12'},...
%     'experiment','vonMisesPrior',...
%     'models',{'model_Bayes_WJM'},'AIC','distanceToPrior=0');

%AIC (log10 scaled)
%%---------------
% 
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10'...
%     'sub11','sub12'},...
%     'experiment','vonMisesPrior',...
%     'dataPath',pathLab,...
%     'models',{'model_CompDiv',...
%     'model_CompDiv_withCard',...
%     'model_Bayes_MAP',...
%     'model_Bayes_Sampling',...
%     'model_Bayes_Sampling_withCard',...
%     'model_Bayes_MAP_FatTailPrior',...
%     'model_Bayes_Sampling_FatTailPriorAndLLH',...
%     'model_Bayes_Sampling_FatTailPrior_ChWT_NoCard',...
%     'model_Bayes_MAP_FatTailPrior_withCard',...
%     'model_Bayes_MAP_withCard'},'AIC','distanceToPrior=90');
% 
%Cross-val R^2 (completed)
%------------------------
% [loglall,AIC,numfitPthisModel,numTrial]=SLdrawAICperCondition({'sub01',...
%     'sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09',...
%     'sub11','sub12'},...
%     'experiment','vonMisesPrior',...
%     'dataPath',pathLab,...
%     'models',{'model_CompDiv',...
%     'model_Bayes_MAP',...
%     'model_Bayes_Sampling'},'CrossValR2');
%
%
%AIC combined vonMises & bimodal priors experiment (completed)
%------------------------------------------------------------
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03'},...
%     'experiment','vonMisesPrior','bimodalPrior',...
%     'dataPath','/Users/steeve/data/dataPsychophy/CombinedExperiments',...
%     'models',{'model_Bayes_MAP','model_Bayes_Sampling'},'AIC');
%
%
%Cross-val maxlogl (completed)
%----------------------------
% [loglall,AIC,numfitPthisModel,numTrial]=SLdrawAICperCondition({'sub01',...
%     'sub02','sub03'},...
%     'experiment','bimodalPrior',...
%     'dataPath','/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%     'models',{'model_Bayes_MAP','model_Bayes_Sampling'},'CrossValmaxlogl');
%
%
%Bimodal Priors experiment (completed)
%-------------------------------------
% [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition({'sub01',...
%     'sub02','sub03'},...
%  'experiment','bimodalPrior',...
%  'dataPath','/Users/steeve/data/dataPsychophy/Exp02_Bimodal_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2',...
%  'models',{'model_CompDiv','model_Bayes_MAP','model_Bayes_Sampling','model_Bayes_MAP_withCard'},'AIC')

function [loglall,AIC,numfitPthisModel,numTrial] = SLdrawAICperCondition(subjects,varargin)

%call for help
if ieNotDefined('subjects')
    help SLdrawAICperCondition
    return
end

%delete varargs which is a common variable across function and may cause
%problems between two functions if called by the two functions.
varargs = varargin;
clear varargin

%check and set data path
if sum(strcmp(varargs,'dataPath'))
    dataPath = varargs{find(strcmp(varargs,'dataPath'))+1};
    
    %check path
    if SLexistPath(dataPath)
    else
        warning('I couldnt find your path...please enter it manually.')
        dataPath = uigetdir;
    end
        cd(dataPath)
else
    fprintf('%s \n','(SLdrawAICperCondition) Please set dataPath')
    dataPath = uigetdir(cd,'Pick up project folder: e.g., data/dataPsychophy/Exp01.../');
    fprintf('%s \n','(SLdrawAICperCondition) Thanks. Your data path was registered.')
    cd(dataPath)
end

%(case von Mises or bimodal prior experiment)
if sum(strcmp(varargs,'experiment'))
    if sum(strcmp(varargs,'vonMisesPrior')) || ...
            sum(strcmp(varargs,'bimodalPrior'))
        myExperiment = varargs(find(strcmp(varargs,'experiment'))+1);
    end
end

%(case combined von Mises and bimodal prior)
if sum(strcmp(varargs,'experiment'))
    if sum(strcmp(varargs,'vonMisesPrior')) && ...
            sum(strcmp(varargs,'bimodalPrior'))
        myExperiment = {'vonMisesPrior','bimodalPrior'};
    end
end

%DeltaI or logDeltaI
if sum(strcmp(varargs,'logDeltaI'))==0 && sum(strcmp(varargs,'DeltaI')) ==0
    deltaOrLogDelta = input('(SLdrawAICperCondition) Do you want to plot DeltaI:1 or logDeltaI:2 ?');
elseif sum(strcmp(varargs,'DeltaI'))==1
    deltaOrLogDelta = 1;
elseif sum(strcmp(varargs,'logDeltaI'))==1
    deltaOrLogDelta = 2; 
end

%(case average deltaI over subjects)
if sum(strcmp(varargs,'subjectAverage'))
    myAnalyses = 'subjectAverage';
else
    myAnalyses = [];
end

%check analyses
if deltaOrLogDelta==1
   fprintf('I will plot the Delta I for each subject and model \n') 
elseif deltaOrLogDelta==2
    fprintf('I will plot the log (deltaI) for each subject and model \n') 
end

%load log likelihood data from modelfit folder and calculate mean log
%likelihood for the experimental condition of interest. Those data are
%obtained from fitting each model to the data
%--------------------------------------------------------------------
%set the directory where you'll find model fit data
%case we are already in data directory go to parent
a = pwd;

%(case AIC from fminsearch fit)
if sum(strcmp(varargs,'AIC'))==1
    
    %go to directory with AIC data
    if strcmp(a(end-4:end),'/modelfit/AIC')
        pathFold = pwd;
    else
        %otherwise go to /modelfit directory
        pathFold = strcat(pwd,'/modelfit/AIC');
    end
    
    %(case AIC from cmaes fit)
elseif sum(strcmp(varargs,'AICcmaes'))==1    
    
    %go to directory with AIC data
    if strcmp(a(end-4:end),'/modelfit/AICcmaes')
        pathFold = pwd;
    else
        %otherwise go to /modelfit directory
        pathFold = strcat(pwd,'/modelfit/AICcmaes');
    end
    
    %(case we want to look at cross-validated R2 to compare models)
    %--------------------------------------------------------------
elseif sum(strcmp(varargs,'CrossValR2'))==1
    
    %go to directory with Cross-validated R^2 data
    if strcmp(a(end-4:end),'/modelfit/CrossValR2/trialData/2sets')
        pathFold = pwd;
    else
        pathFold = strcat(pwd,'/modelfit/CrossValR2/trialData/2sets');
    end
    
    %(case we want to look at cross-validated maxlogl to compare models)
    %-------------------------------------------------------------------
elseif sum(strcmp(varargs,'CrossValmaxlogl'))==1
    
    %go to directory with Cross-validated R^2 data
    if strcmp(a(end-4:end),'/modelfit/CrossValmaxlogl/trialData/2sets')
        pathFold = pwd;
    else
        pathFold = strcat(pwd,'/modelfit/CrossValmaxlogl/trialData/2sets');
    end
end

%status
fprintf('%s \n','(SLdrawAICperCondition) I am getting AIcs in this folder: ',pathFold)

%get the names and informations about model folders found in the directory
d = dir([pathFold '/*model*']);

%models
if sum(strcmp(varargs,'models'))
    models = varargs{find(strcmp(varargs,'models'))+1};
else
    fprintf('%s \n','Please enter the name of two models...')
end

%find the subfolders to analyse
if size({d.name},2) < numel(models)
    
    %display wrong models
    for i=1:numel(models)
        for j=1:size({d.name},2)
            wrongModels{i,j} = intersect(d(j).name,models(i));
        end
    end
    pos = sum(SLisCellEmpty(wrongModels),2)==size({d.name},2);
    wrongModels = {models{pos}};
    numWrong = num2str(size(wrongModels,2));
           
    %warning
    fprintf('%s \n',[numWrong,' models you entered do not exist.'])
        
    for i = 1 : size(wrongModels,2)
        fprintf('\n %s \n',['  - ',wrongModels{i}])
    end
    
    fprintf('\n %s','Enter at least 2 of those models:',...
        '- model_Bayes_MAP',...
        '- model_Bayes_Sampling',...
        '- model_Bayes_Sampling_withCard',...
        '- model_Bayes_MAP_FatTailPrior',...
        '- model_Bayes_Sampling_FatTailPriorAndLLH',...
        '- model_Bayes_Sampling_FatTailPrior_ChWT_NoCard',...
        '- model_Bayes_MAP_FatTailPrior_withCard',...
        '- model_Bayes_MAP_withCard',...
        '- model_CompDiv',...
        '- model_CompDiv_withCard',...
        '- model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails')
    fprintf('\n %s','')

    keyboard
end
for ii = 1 : numel(models)
    imodel(ii) = find(strcmp({d.name},models(ii)));
end

%subjects
subjectsAIC = subjects;

%get the names of the subfolders
Folds.name = {d(imodel).name}';

%count the number of subfolders
Folds.nb = numel(Folds.name);

%initialize variables
logL_pertrialBkp = [];
FA1 = [];
FA2 = [];
FA3 = [];
TheSub = [];
TheModels = [];
numfitP = [];
loglall         =[];
AIC             =[];
numfitPthisModel=[];
numTrial        =[];


%case von Mises or bimodal prior
%-------------------------------
if sum(strcmp(myExperiment,'vonMisesPrior')) && sum(strcmp(myExperiment,'bimodalPrior'))==0 ...
        || sum(strcmp(myExperiment,'vonMisesPrior'))==0 && sum(strcmp(myExperiment,'bimodalPrior'))==1
    
    %---
    %AIC
    %---
    if sum(strcmp(varargs,'AIC')) || sum(strcmp(varargs,'AICcmaes'))
        
        %Create a data matrix of data and variables of interest
        %loop over the model folders
        for j = 1 : Folds.nb
            
            %clear up variables
            clear nb_of_fitp
            clear fitP
            clear o

            %set this model folder
            thisModelPath = strcat([pathFold,'/',Folds.name{j}]);
            
            %get the names and informations about subjects folders found in the
            %directory   
            %ignore missing subjects
            for i = 1 : length(subjectsAIC)
                ourdir = dir([thisModelPath '/*' num2str(subjectsAIC{i}) '*']);
                if ~isempty(ourdir)
                    notmissing(i) = 1;
                else
                    notmissing(i) = 0;
                end
            end
            pos = find(notmissing==1);
            clear dsub
            for i = 1 : length(subjectsAIC)
                dsub{i} = dir([thisModelPath '/*' num2str(subjectsAIC{i}) '*']); 
            end
            
            %get subjects
            dsub = [dsub{:}];
            
            %get the number of subjects
            numSub = length(dsub);
            
            %get the names of subjects' subfolders
            SubFolds.name = {dsub(1:1:numSub).name}';
            
            %display models and subjects
            fprintf('%s \n', ['Model is: ',models{j}])
                        
            %collect data for each subject
            for k = 1 : numSub
                
                %get and store subject
                subject = SubFolds.name{k};
                fprintf('%s \n',subject)

                %set this subject path and go there
                thisModelsandSubPath = [thisModelPath,'/',subject];
                cd(thisModelsandSubPath)
                
                %load and list data
                datalisting = dir('*datafit*.mat');
                
                %warning
                if isempty(datalisting)
                    sprintf(['(SLdrawAICperCondition) No data for this subject yet. You should delete his folder'])
                    keyboard
                end

                %load data
                load(datalisting.name);
                
                %number of trials
                %check that logL per trials is in the working space (ws)
                %make sure the variable name is correct
                c0 = who('Logl_pertrialBestfit');
                if ~isempty(c0)
                   logL_pertrialBestfit = Logl_pertrialBestfit;
                   clear Logl_pertrialBestfit
                end
                
                %check that the correct variable is in the ws
                c = who('logL_pertrialBestfit');
                if isempty(c)
                    if isempty(who('o'))
                      fprintf('%s \n','(SLdrawAICperCondition) "logL_pertrialBestfit" was not found in the working space')
                      keyboard
                    else
                        %case only sum of log likelihood is found
                        if numel(o.logl)==1
                            logL_pertrialBestfit = nan(length(o.posC),1);    
                            logL_pertrialBestfit(1) = o.logl; 
                            fprintf('%s \n','(SLdrawAICperCondition) Sum logl only was found for this subject. Analyses with loglpertrial will not be run.') 
                        end    
                    end
                end

                %get and store log likelihood of each trial data given the model
                %logL_pertrial = [logL_pertrial ; Logl_pertrialBestfit];
                logL_pertrialBkp = [logL_pertrialBkp ; logL_pertrialBestfit];
               
                %numTrials = length(output.Logl_pertrialBestfit);
                numTrials = length(logL_pertrialBestfit);
                clear logL_pertrialBestfit;
                clear Logl_pertrialBestfit;

                %models and subjects' name
                TheModels = [TheModels ; repmat({models{j}},numTrials,1)];
                TheSub = [TheSub; repmat({subject},numTrials,1)];
                
                %number of model parameters
                %case 'fitP' exists
                if ~isempty(who('fitP'))
                    nb_of_fitp = sum(~isnan(fitP.p));
                elseif ~isempty(who('o'))
                    nb_of_fitp = sum(~isnan(o.fitP));
                else
                    fprintf('%s \n','(SLdrawAICperCondition) Could not find the "fitp"')
                     keyboard
                end
                numfitP  = [numfitP; repmat({nb_of_fitp},numTrials,1)];
                
                %get and store experimental conditions (e.g., motion coherence, prior strength, etc.,...)
                %get factor 1
                %get factor 2
                %(case von Mises prior)
                %---------------------
                c = who('output');
                
                %prior strength
                if strcmp(myExperiment,'vonMisesPrior') && strcmp(myExperiment,'bimodalPrior')==0
                                 
                    %check if the variable output exists
                    if isempty(who('output')) && isempty(who('o'))

                         output.StimStrengthVM = StimStrength;
                         output.pstdVM         = pstd;
                         output.dispVM         = disp(:);

                    else

                        %check if the fields for task conditions: coh, pstd or
                        %disp exist
                        if ~isempty(who('output')) && isfield(output,'coh') && isfield(output,'pstd')  && isfield(output,'disp')

                            output.StimStrengthVM = output.coh;
                            output.pstdVM = output.pstd;
                            output.dispVM = output.disp(:);
                        
                        %case output exists and the factors                                
                        elseif ~isempty(who('output')) && isfield(output,'StimStrength') && isfield(output,'pstd')  && isfield(output,'disp')
                            
                            output.StimStrengthVM = output.StimStrength;
                            output.pstdVM = output.pstd;
                            output.dispVM = output.disp(:);
                            
                        elseif ~isempty(who('output')) && isfield(output,'stimStrength') && isfield(output,'pstd')  && isfield(output,'disp')
                            
                            output.StimStrengthVM = output.stimStrength;
                            output.pstdVM = output.pstd;
                            output.dispVM = output.disp(:);
                            
                        %case o exist but not the factors    
                        elseif ~isempty(who('o')) && ~isfield(o,'coh') && ~isfield(o,'pstd')  && ~isfield(o,'disp')
                            
                            output.StimStrengthVM = o.uniqCond(o.posC,2);
                            output.pstdVM         = o.uniqCond(o.posC,1);
                            output.dispVM         = o.uniqCond(o.posC,3);

                        end
                    end
                    
                    %update the task conditions over subjects and models
                    FA1 = [FA1 ; output.StimStrengthVM];
                    FA2 = [FA2 ; output.pstdVM];
                    FA3 = [FA3 ; output.dispVM(:)];
                    FA1nm = 'coh';
                    FA2nm = 'Pstd';
                    FA3nm = 'disp';
                end
                
                %(case bimodal prior)
                %--------------------
                %the second factor is the priors' modes
                if strcmp(myExperiment,'bimodalPrior') && strcmp(myExperiment,'vonMisesPrior')==0
                    
                    %check if output variable exists
                    if isempty(c)==1
                         output.StimStrengthBim = StimStrength;
                         output.pstdBim = pstd;
                         output.priorModesBim = priorModes;
                         output.dispBim = disp(:) ;
                    end
                    
                    FA1 = [FA1 ; output.StimStrengthBim];
                    FA2 = [FA2 ; output.priorModesBim];
                    FA3 = [FA3 ; output.dispBim(:)];
                    FA1nm = 'coh';
                    FA2nm = 'priormodes';
                    FA3nm = 'disp';
                end
                
                clear output
            end
        end
        
        %clear databank found
        databank.data = [];
        
        %create our databank
        %enter the names of the variable in the columns
        databank.fitnm=[
            {'logL_pertrial'},...
            {'coh'},...
            {FA2nm},...
            {'sample_dir'},...
            {'subjects'},...
            {'TheModels'},...
            {'numfitP'}];
        
        for i = 1 : size(FA2,1)
            FA2bkp{i} = FA2(i,:);
        end
        
        %store the data in each column
        databank.fitdata = [
            num2cell(logL_pertrialBkp) ...
            num2cell(FA1) ...
            FA2bkp' ...
            num2cell(FA3) ...
            TheSub ...
            TheModels ...
            numfitP];
        
        %warning
        if length(databank.fitnm) ~= size(databank.fitdata,2)
            fprintf('%s \n','There is a problem in the databank')
            keyboard
        end
        
        %backup the databank
        save('datbank_fit','databank');
        
        %data are saved in a file called datbank_fit in the parent directory
        cd(pathFold)
        
        %AICs for conditions of interest (that maximize model comparison)
        %----------------------------------------------------------------
        %status
        fprintf('%s \n','Calculating AIC...')
        
        %(case von Mises prior experiment)
        %---------------------------------
        %check that experiment is input
        if sum(strcmp(varargs,'experiment'))
            if strcmp(varargs(find(strcmp(varargs,'experiment'))+1),...
                    'vonMisesPrior')
                
                %status
                fprintf('%s \n','(SLdrawAICperCondition) Experiment: von Mises prior...')
                
                %experimental conditions of interest (where the model
                %predictions diverge most)
                %--------------------------------------------------------
                DispDir = cell2mat(databank.fitdata(:,(strcmp(databank.fitnm,...
                    'sample_dir'))==1));
                
                %calculate motion direction distance to the prior
                Dist2Prior = SLvectors2signedAngle(DispDir,225,'polar');
                
                %select data for all motion directions
                %myCondition = abs(Dist2Prior)>=0;
                
                %select data when motion direction was at least 
                %'distanceToPrior' deg distant from the prior.
                if sum(strcmp(varargs,'distanceToPrior'))==1
                    
                    pos = find(strcmp(varargs,'distanceToPrior')==1);
                    myDis2Prior= varargs{pos+1};
                    myCondition = abs(Dist2Prior)>=myDis2Prior;
                    
                    %status
                    fprintf(['(SLdrawAICperCondition) I am selecting data where displayed direction is at least ',num2str(myDis2Prior),...
                        ' distant from the prior to calculate the AIC\n'])
                else
                    
                    %get entire dataset
                    myCondition = abs(Dist2Prior)>=0;
                    
                    %status
                    fprintf(['(SLdrawAICperCondition) I am selecting the entire dataset to calculate the AICs'])
                    
                end
                databank.fitdatamyCondition = databank.fitdata(myCondition,:);
                save('datbank_fit','databank');
                
                %AICs
                %----
                %AIC is calculated like this "AIC=2k - 2ln(L)". We calculate a
                %matrix of AIC values for each subject (rows) and each models
                %(column).
                %ln(L) is the log likelihood of the entire dataset given each
                %model. It is the sum of the log likelihood of each trial-data.
                %ln(L) = ln(L(data1/model1)) + ln(L(data2/model1)) + ....
                
                %logL of data given models
                logL_pertrialMyC = cell2mat(databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'logL_pertrial'))==1));
                
                %models' trials
                ModelThisTrial = databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'TheModels'))==1);
                
                %models
                thisModel = unique(ModelThisTrial,'stable');
                numModels = length(thisModel);
                
                %subjects' trials
                SubThisTrial = databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'subjects'))==1);
                
                %number of models
                numfitP = cell2mat(databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'numfitP'))==1));
                
                %last subject
                uniqSubj = unique(SubThisTrial);
                numTotalSub = numel(uniqSubj);
                
                %this model's trials
                %get number of fit parameter
                %get subjects tested for this model
                %SubThisModel=nan(20,numModels);
                AIC = nan(length(numModels),length(numTotalSub));
                loglall = nan(length(numModels),length(numTotalSub));
                numTrial = nan(length(numModels),length(numTotalSub));
                numfitPthisModel = nan(numModels,1);
                
                for i = 1 : numModels
                    
                    %model trials
                    ThisModelTrials = ismember(ModelThisTrial,thisModel{i});
                    
                    %nb of fitP
                    nfitP = unique(numfitP(ThisModelTrials),'stable');
                    
                    %warning
                    if numel(nfitP)~=1
                        sprintf(['(SLdrawAICperCondition) ""',thisModel{i},'" cannot',...
                            ' be assigned more than 1 number of fit parameters'])
                        keyboard
                    end
                    numfitPthisModel(i) = unique(numfitP(ThisModelTrials),'stable');
                    SubThisModel = unique(SubThisTrial(ThisModelTrials));
                    
                    %AICs
                    %loop over all subjects
                    for j = 1 : numTotalSub
                        if sum(strcmp(SubThisModel,uniqSubj{j}))==1
                            
                            %this subject's trials
                            ThisSubTrials = ismember(SubThisTrial,uniqSubj{j});
                            
                            %AIC and DeltaI

                            if ~isnan(sum(logL_pertrialMyC(ThisModelTrials & ThisSubTrials)))

                                %case logl per trial was found in the dataset and no missing
                                loglall(i,j) = sum(logL_pertrialMyC(ThisModelTrials & ...
                                    ThisSubTrials));

                            else

                                %case only sumlogl was found in dataset (not logl per trial)
                                loglall(i,j) = nansum(logL_pertrialMyC(ThisModelTrials & ...
                                    ThisSubTrials));
                                fprintf('%s \n','(SLdrawAICperCondition) only sumlogl (not logl per trial) was found in dataset and was used then to calculate AIC.')

                            end

                            AIC(i,j) = 2*numfitPthisModel(i)-2*loglall(i,j);
                            
                            %number of trials eacn condition
                            numTrial(i,j) = sum(ThisModelTrials & ThisSubTrials);
                        else
                            %missing data
                            loglall(i,j) = NaN;
                            AIC(i,j) = NaN;
                        end
                    end
                end
                
                %models and subjects names
                modelnm = unique(ModelThisTrial,'stable');
                
                %draw summary Table
                rowLabels = unique(TheSub);
                
                %col
                columnLabels = ['"AIC (2*(nP -logl)"' ; modelnm ; ...
                    '"nb_of_fitP"' ; modelnm ;...
                    '"logl"';modelnm];
                
                %data
                rowTitles = nan(1,length(rowLabels));
                dataTable = [rowTitles; AIC; ...
                    rowTitles; numfitPthisModel(:,ones(1,numTotalSub));...
                    rowTitles; loglall];
                
                %table
                fprintf('(SLdrawAICperCondition) I am printing the data table \n')
                SLprintTable(dataTable,columnLabels,rowLabels)
                
                %draw AICs (model,subjects)
                %figure('color','w')
                
                %color bars and text
                a = 0.9./numModels;
                colors = repmat([a:a:1]',1,3);               
                txtc = 1 - repmat([a:a:1]',1,3);
                
                %position bars on x axis
                for i = 1 : numModels
                    spaceBtwSub = 5;
                    x1(i,:) = [1:numModels+spaceBtwSub:(numModels+spaceBtwSub)*numTotalSub - numModels]+(i-1);
                end
                
                %draw AIC for subject and models
%                 fprintf('(SLdrawAICperCondition) I am drawing AICs for each subject and models \n')
%                 
%                 for i = 1 : numTotalSub
%                     for j = 1 : numModels
%                         hold all
% 
%                         %bars (Stocker and Simoncelli, 2006,NN)
%                         b(j) = bar(x1(j,i),AIC(j,i),...
%                             'edgecolor','none',...
%                             'facecolor',colors(j,:),...
%                             'displayname',modelnm{j});
%                         
%                         %annotate
%                         %text(x1(j,i),1.07*AIC(j,i),num2str(fix(AIC(j,i)*100)/100),...
%                         %   'HorizontalAlignment','center','fontsize',10)
%                     end
%                 end
%                 %ylim([0 1.05*max(AIC(:))])
%                 xlim([0 x1(end)+1])
%                 ylabel('AIC','fontsize',14)
%                 drawPublishAxis
                
                %draw DeltaI (AIC - min AIC)
                %---------------------------
                fprintf('(SLdrawAICperCondition) I am drawing Delta Is for each subject and models \n')
                
                %Be careful of DeltaI close to 0 it is
                %difficult to interprete.
                figure('color','w','position',[469 309 1012 241])
                hold all
                maxDeltaI = nan(numTotalSub,1);
                SubID     = nan(numTotalSub,1);
                deltaIBKP = nan(numModels,numTotalSub);
                
                %case sorting subjects by best model
                if sum(strcmp(varargs,'sortByBestModel'))==1
                    if numModels == 2
                        
                        fprintf('(SLdrawAICperCondition) Sorting subjects by best model \n')
                        [i,pos] = sort(AIC(1,:) - AIC(2,:));
                        AIC = AIC(:,pos);
                        
                    end
                else 
                    pos = 1:numTotalSub;
                end
                
                %get sub ID
                for i = 1 : numTotalSub 
                    SubID(i) = str2double(uniqSubj{i}(4:end));    
                end
                SubID = SubID(pos);
                
                %Plot delta I or log Delta I for each subject and model
                for i = 1 : numTotalSub
                    
                    %suject ID
                    %SubID(i) = str2double(uniqSubj{i}(4:end));
                    
                    %models
                    for j = 1 : numModels
                        
                        %ax
                        %ca(i)=subplot(1,numTotalSub,i);
                        hold all
                        
                        %when data
                        if ~isnan(AIC(j,i))
                            
                            deltaI = AIC(j,i) - min(AIC(:,i));
                            
                            %backup
                            deltaIBKP(j,i) = deltaI;

                            %logDeltaI
                            if deltaOrLogDelta == 2
                                %bars (Stocker and Simoncelli, 2006,NN)

                                %add 1 just for visbility.
                                space = 0.03;
                                
                                %plot
                                %trick because delta = 0 (best model) is -inf
                                %if logged. Add space to visualize delta = 0.
                                deltaI(deltaI==0) = 1;
                                b(j) = bar(x1(j,i),log10(deltaI) + space,...
                                    0.95,...
                                    'edgecolor','none',...
                                    'facecolor', colors(j,:),...
                                    'displayname', modelnm{j}(7:end));
                                
                                %annotate
                                text(x1(j,i),(log10(deltaI)+space)./2,num2str(SLround(log10(deltaI),0.1)),...
                                    'color',txtc(j,:),'verticalAlignment','middle',...
                                    'rotation',90,'fontsize',12)
                                
                                %plot deltaI
                            elseif deltaOrLogDelta == 1 
                                
                                %visibility
                                space = 10;
                                
                                b(j) = bar(x1(j,i),deltaI+space,...
                                    0.95,...
                                    'edgecolor','none',...
                                    'facecolor', colors(j,:),...
                                    'displayname', modelnm{j}(7:end));
                                
                                %annotate
                                text(x1(j,i),(deltaI+space)./2,num2str(SLround(deltaI,0.1)),...
                                    'color',txtc(j,:),'verticalAlignment','middle',...
                                    'rotation',90,'fontsize',12)                                                               
                            else
                                fprintf(['(SLdrawAICperCondition) ' ...
                                         'Please set "DeltaI" or ' ...
                                         '"logDeltaI" as an argument'])
                                keyboard
                            end
                        end
                        
                        %graphics
                        maxDeltaI(i) = max(AIC(:,i) - min(AIC(:,i)));
                        
                    end
                    
                    %set labels starting 0 (with 10 added for visibility)
                    set(gca,'ytick',10:100:max(maxDeltaI),'yticklabel',0:100:max(maxDeltaI))
                    
                    %graphics
                    %xlim([0 numModels+1])
                    %set(gca,'fontsize',11)
                    %set(gca,'xticklabel',[])
                    %if i~=1
                    %    set(gca,'ycolor','w')
                    %end
                    
                    %labels
                    ylabel('$\Delta$ I','fontsize',14,'interpreter','latex')
                    xlabel('Subject','fontsize',14)
                end
         
                %label subjects
                set(gca,'xtick',x1(1,:),'xticklabel',SubID)
                xlim([0 max(x1(:))+1])
                
                %indicate AIC difference of 2 (significant evidence in
                %favor to best model)
                %visibility
                space = 10;                
                plot(linspace(1,max(x1(:))),linspace(log10(2)+space,log10(2)+space),...
                   '--k')

                %SLConventionUp
                %legend(b,'Location','EastOutside')
                
                
                %(case subject average delta I)
                %...............................
                if sum(strcmp(myAnalyses,'subjectAverage'))
                    
                    %arrange data
                    %data = deltaIBKP(:);
                    %x1 = repmat(1:numModels,1,numTotalSub)';
                    %x2 = repmat(1:numTotalSub,numModels,1);
                    %x2 = x2(:);
                    
                    %stats
                    %myStats = SLmakeStat(data,x1,x2);
                    
                    %mean and variability over sub
                    meandeltaI = nanmean(deltaIBKP,2);
                    if size(deltaIBKP,2)>1
                        stddeltaI = nanstd(deltaIBKP');
                        semdeltaI = sem(deltaIBKP,1);
                    else
                        stddeltaI = zeros(length(deltaIBKP),1);
                        semdeltaI = zeros(length(deltaIBKP),1);
                    end
                    
                    %sort
                    [sort_meandeltaI,pos] = sort(meandeltaI);
                    sort_stddeltaI = stddeltaI(pos);
                    sort_semdeltaI = semdeltaI(pos);
                    sort_models = models(pos);
                    sort_modelnm = modelnm(pos);
                    
                    %draw mean deltaI over sub
                    figure('color','w');
                    for i =  1 : numModels
                        hold all
                        
                        %bar
                        sp(i) = bar(i,sort_meandeltaI(i),...
                            0.95,...
                            'edgecolor','none',...
                            'facecolor', colors(i,:),...
                            'displayname', sort_modelnm{i}(7:end));
                        
                        %errorbar
                        SLerrorbar(i,sort_meandeltaI(i),'yError',...
                            sort_semdeltaI,'Symbol=-')
                        
                        %values 
                        text(i,sort_meandeltaI(i)./4,num2str(SLround(sort_meandeltaI(i),0.1)),...
                                'color',txtc(i,:),'verticalAlignment','middle',...
                                'rotation',90,'fontsize',12)
           
                    end
                    
                    %indicate AIC difference of 2 (significant evidence in
                    %favor to best model)
                    hold all
                    %plot(linspace(1,numModels),linspace(2 + sort_meandeltaI(1),2 + sort_meandeltaI(1)),...
                    %    '--r')
                    %labels
                    ylabel('Subjects average $\Delta$ I','fontsize',14,'interpreter','latex')
                    xlabel('Models','fontsize',14)
                    
                    axis tight
%                     SLConventionUp
                    legend(sp,'Location','EastOutside')

                end
            end
        end
        
        
        %-------------------------------
        %(case bimodal prior experiment)
        %-------------------------------
        %check that experiment is input
        if sum(strcmp(varargs,'experiment'))
            if strcmp(varargs(find(strcmp(varargs,'experiment'))+1),...
                    'bimodalPrior')
                
                %status
                fprintf('%s \n','Experiment: Bimodal prior...')
                
                %set experimental conditions of interest (where the model
                %predictions diverge most)
                %--------------------------------------------------------
                DispDir = cell2mat(databank.fitdata(:,...
                    (strcmp(databank.fitnm,'sample_dir'))==1));
                
                %calculate distance of displayed direction to the prior
                %counterclock mode at each trial
                %when the distance between the first and second mode is <0 the
                %first mode is more counterclockwise than the second and when it is
                %>0, the second mode is more clockwise than the first. Now get the
                %coordinates of the counterclock mode (left)
                %Dist2Prior=vectors2signedAngle(DispDir,225);
                distModes = SLvectors2signedAngle(FA2(:,1),FA2(:,2),'polar');
                leftMode(distModes<0) = FA2(distModes<0,1);
                leftMode(distModes>0) = FA2(distModes>0,2);
                rightMode(distModes<0) = FA2(distModes<0,2);
                rightMode(distModes>0) = FA2(distModes>0,1);
                Dist2leftMode = SLvectors2signedAngle(DispDir,leftMode,'polar');
                Dist2rightMode = SLvectors2signedAngle(DispDir,rightMode,'polar');
                
                %select direction that are not in between the two modes where
                %model predictions differ most (3 modes for competition vs.
                %modes for Bayesian). The furthest the best (the tails are 30?
                %away from the modes)
                myCondition = Dist2leftMode<=0 | Dist2rightMode>=0;
                
                %all data
                %myCondition=abs(distModes)>=0;
                
                %save data
                databank.fitdatamyCondition=databank.fitdata(myCondition,:);
                save('datbank_fit','databank');
                
                %AICs
                %----
                %AIC is calculated like this "AIC = 2k - 2ln(L)". We calculate a
                %matrix of AIC values for each subject (rows) and each models
                %(column).
                %note that we use expected loglikelihood over trials instead of
                %likelihood that we couldn't evaluate because the factorial
                %product of so many probabilities gives 0 as a likelihood.
                %AICc is better to use but because it requires likelihood, we
                %cannot use it.
                
                %logL data
                logL_pertrialMyC=cell2mat(databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'logL_pertrial'))==1));
                
                %models' trials
                ModelThisTrial = databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'TheModels'))==1);
                
                %models
                thisModel = unique(ModelThisTrial,'stable');
                numModels = length(thisModel);
                
                %subjects' trials
                SubThisTrial = databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'subjects'))==1);
                
                %number of models
                numfitP = cell2mat(databank.fitdatamyCondition(:,...
                    (strcmp(databank.fitnm,'numfitP'))==1));
                
                %last subject
                uniqSubj=unique(SubThisTrial);
                numTotalSub=numel(uniqSubj);
                
                %this model's trials
                %number of fit parameter
                %subjects tested for this model
                %SubThisModel=nan(20,numModels);
                AIC = nan(length(numModels),length(numTotalSub));
                loglall = nan(length(numModels),length(numTotalSub));
                numTrial = nan(length(numModels),length(numTotalSub));
                for i = 1 : numModels
                    ThisModelTrials = ismember(ModelThisTrial,thisModel{i});
                    numfitPthisModel(i) = unique(numfitP(ThisModelTrials),'stable');
                    SubThisModel = unique(SubThisTrial(ThisModelTrials));
                    
                    %AICs
                    for j = 1 : numTotalSub
                        if sum(strcmp(SubThisModel,uniqSubj{j}))==1
                            
                            %this subject's trials
                            ThisSubTrials = ismember(SubThisTrial,uniqSubj{j});
                            
                            %mean log likelihood and AIC
                            loglall(i,j) = sum(logL_pertrialMyC(ThisModelTrials & ...
                                ThisSubTrials));
                            AIC(i,j) = 2*numfitPthisModel(i)-2*loglall(i,j);
                            
                            %number of trials
                            numTrial(i,j) = sum(ThisModelTrials & ThisSubTrials);
                        else
                            %missing data
                            loglall(i,j) = NaN;
                            AIC(i,j) = NaN;
                        end
                    end
                end
                
                %draw AICs for the different models for each subject
                %---------------------------------------------------
                figure('color','w')
                colors=[0 0 0;
                    .5 .5 .5;
                    .7 .7 .7;
                    .8 .8 .8;
                    .9 .9 .9];
                
                %position data
                for i = 1:numModels
                    x1(i,:)=[1:numModels+1:(numModels+1)*numTotalSub-numModels]+(i-1);
                end
                
                %get models and subjects names
                modelnm = unique(ModelThisTrial,'stable');
                for i = 1 : numTotalSub
                    for j = 1 : numModels
                        hold all
                        b(j) = bar(x1(j,i),AIC(j,i),...
                            'edgecolor','none',...
                            'facecolor',colors(j,:),...
                            'displayname',modelnm{j});
                        
                        %annotate
                        text(x1(j,i),1.07*AIC(j,i),num2str(fix(AIC(j,i)*100)/100),...
                            'HorizontalAlignment','center','fontsize',10)
                    end
                end
                ylim([0 1.01*max(AIC(:))])
                xlim([0 x1(end)+1])
                ylabel('AIC','fontsize',14)
                drawPublishAxis
                
                %draw DeltaI (AIC - min AIC)
                figure('color','w')
                for i = 1 : numTotalSub
                    SubID(i) = str2double(uniqSubj{i}(4:end));
                    
                    for j = 1 : numModels
                        hold all
                        
                        %if data
                        if ~isnan(AIC(j,i))
                            
                            %bars (Stocker and Simoncelli, 2006,NN)
                            %add 1 just for visbility.
                            %same axes
                            deltaI = AIC(j,i) - min(AIC(:,i));
                            space = 0.03;
                            
                            %plot
                            %trick because delta = 0 (best model) is -inf
                            %if logged. Add space to visualize delta = 0.
                            deltaI(deltaI==0) = 1;
                            b(j) = bar(x1(j,i),log10(deltaI) + space,...
                                0.95,...
                                'edgecolor','none',...
                                'facecolor', colors(j,:),...
                                'displayname', modelnm{j}(7:end));
                        else
                            %if missing data
                            scatter(j,0,100,...
                                'k');
                        end
                        
                        %graphics
                        maxDeltaI(i) = max(AIC(:,i)-min(AIC(:,i)));
                    end
                    %labels
                    ylabel('Log $\delta$ I','fontsize',14,'interpreter','latex')
                    xlabel('Subject','fontsize',14)
                end
                legend(regexprep(b,'_',' '))
                
                %label subjects
                set(gca,'xtick',x1(1,:),'xticklabel',SubID)
                xlim([0 max(x1(:))+1])
                
                %indicate AIC difference of 2. (significance evidence for the
                %best model)
                plot(linspace(1,max(x1(:))),linspace(log10(2)+space,log10(2)+space),...
                    '--k')
                %SLConventionUp(gcf)
            end
        end
    end
    
    
    
    %--------------------------------------
    %case log likelihood for all conditions
    %--------------------------------------
    %subjects, models, and experimental factors
    if sum(strcmp(varargs,'loglEachCondition'))==1
        
        %logL data
        logL_pertrialMyC=cell2mat(databank.fitdata(:,(strcmp(databank.fitnm,...
            'logL_pertrial'))==1));
        
        %experimental factors
        FA1thisTrial=cell2mat(databank.fitdata(:,(strcmp(databank.fitnm,...
            'coh'))==1));
        FA2thisTrial=cell2mat(databank.fitdata(:,(strcmp(databank.fitnm,...
            'Pstd'))==1));
        FA3thisTrial=cell2mat(databank.fitdata(:,(strcmp(databank.fitnm,...
            'sample_dir'))==1));
        
        %models' trials
        ModelThisTrial=databank.fitdata(:,(strcmp(databank.fitnm,...
            'TheModels'))==1);
        thisModel=unique(ModelThisTrial,'stable');
        numModels=length(thisModel);
        
        %subjects trials
        SubThisTrial=databank.fitdata(:,(strcmp(databank.fitnm,...
            'subjects'))==1);
        
        %AIC for models and subject
        for i=1:numModels
            
            %model's trials
            ThisModelTrials=ismember(ModelThisTrial,thisModel{i});
            
            %subjects and number of subjects for this model
            SubThisModel=unique(SubThisTrial(ThisModelTrials));
            numSubThisModel=length(SubThisModel);
            
            for j=1:numSubThisModel
                
                %this subject's trials
                ThisSubTrials=ismember(SubThisTrial,...
                    SubThisModel{j});
                
                %mean log likelihood (comparison between models
                %are valid if both models have the same number of parameter)
                logL_thisCondition=logL_pertrialMyC(ThisModelTrials & ...
                    ThisSubTrials);
                
                %sort experimental conditions for this model and this subject
                FA1_thisCondition=FA1thisTrial(ThisModelTrials & ...
                    ThisSubTrials);
                FA2_thisCondition=FA2thisTrial(ThisModelTrials & ...
                    ThisSubTrials);
                FA3_thisCondition=FA3thisTrial(ThisModelTrials & ...
                    ThisSubTrials);
                
                %draw mean log likelihood for experimental each condition
                drawStat(logL_thisCondition,FA3_thisCondition',...
                    FA1_thisCondition',FA2_thisCondition');
                
                %graphics
                xlim([0 360])
                ylim([-5 0])
                title([thisModel{i},' ',SubThisModel{j}],'fontsize',14)
            end
        end
    end
    
    
    
    %-----------------------------------------------
    %(case cross-validated maxlogl to compare models
    %-----------------------------------------------
    if sum(strcmp(varargs,'CrossValR2'))==1
        
        %Create a  matrix of data and variables of interest
        %loop over the model folders
        numSub = numel(subjects);
        numModels = numel(models);
        mysubjects = subjects;
        
        %initialize variable
        meanCrossValR2 = nan(Folds.nb,numSub);
        
        for j=1:Folds.nb
            
            %set this model folder
            thisModelPath=strcat([pathFold,'/',Folds.name{j}]);
            
            %get the names and informations about subjects folders found in the
            %directory
            %dsub=dir([thisModelPath '/*sub*']);
            
            %get the names and informations about input subjects
            %dsub=nan(numSub,1);
            clear dsub
            for i=1:numSub
                dsub(i)=dir([thisModelPath '/*' num2str(mysubjects{i}) '*']);
            end
            
            %get the number of subjects
            numSub=length(dsub);
            
            %get the names of subjects' subfolders
            SubFolds.name={dsub(1:1:numSub).name}';
            
            %display models and subjects
            numModels=numel(models);
            fprintf('%s \n',['Model is: ',models{j}])
            fprintf('%s \n', 'Subjects are:')
            display(SubFolds.name)
            
            %collect data for each subject
            for k=1:numSub
                
                %get and store subject
                subject=SubFolds.name{k};
                
                %set this subject path and go there
                thisModelsandSubPath = [thisModelPath,'/',subject];
                cd(thisModelsandSubPath)
                
                %load and list data
                datalisting = dir('*datafit*.mat');
                load(datalisting.name);
                
                %mean Cross-valudated R-squared
                meanCrossValR2(j,k) = output.CrossVmeanR2;
            end
        end
        
        
        %position data
        for i=1:numModels
            x1(i,:)=[1:numModels+1:(numModels+1)*numSub-numModels]+(i-1);
        end
        
        %colors
        colors=[0 0 0;
            .5 .5 .5;
            .7 .7 .7;
            .8 .8 .8;
            .9 .9 .9];
        
        %status
        fprintf('%s \n','Plotting Cross-validated R2...')
        
        %draw cross validated R-squared
        figure('color','w')
        maxCVR2=nan(numSub,1);
        for i=1:numSub
            b=nan(numModels,1);
            for j=1:numModels
                
                %ax
                hold all
                
                %when there is data
                if ~isnan(meanCrossValR2(j,i))
                    
                    %bars (Stocker and Simoncelli, 2006,NN)
                    b(j)=bar(x1(j,i),meanCrossValR2(j,i),...
                        0.95,...
                        'edgecolor','none',...
                        'facecolor',colors(j,:),...
                        'displayname',models{j}(7:end));
                else
                    %if missing data
                    scatter(j,0,100,...
                        'k');
                end
                
                %graphics
                maxCVR2(i)=max(meanCrossValR2(:));
            end
            
            %graphics
            %labels
            ylabel('Mean cross-validated R-squared','fontsize',14)
            xlabel('Subject','fontsize',14)
        end
        legend(b);
        
        %label subjects
        set(gca,'ytick',(0:0.1:max(maxCVR2)),'yticklabel',...
            (0:0.1:max(maxCVR2)))
        set(gca,'xtick',x1(1,:),'xticklabel',mysubjects)
        xlim([0 max(x1(:))+1])
    end
    
    %------------------------------------------
    %(case cross-validated R2 to compare models
    %------------------------------------------
    if sum(strcmp(varargs,'CrossValmaxlogl'))==1
        
        %Create a  matrix of data and variables of interest
        %loop over the model folders
        numSub = numel(subjects);
        numModels = numel(models);
        mysubjects = subjects;
        
        %initialize variable
        meanCrossValmaxlogl = nan(Folds.nb,numSub);
        
        for j = 1 : Folds.nb
            
            %set this model folder
            thisModelPath = strcat([pathFold,'/',Folds.name{j}]);
            
            %get the names and informations about subjects folders found in the
            %directory
            %dsub=dir([thisModelPath '/*sub*']);
            
            %get the names and informations about input subjects
            %dsub=nan(numSub,1);
            clear dsub
            for i = 1 : numSub
                dsub(i) = dir([thisModelPath '/*' num2str(mysubjects{i}) '*']);
            end
            
            %get the number of subjects
            numSub = length(dsub);
            
            %get the names of subjects' subfolders
            SubFolds.name = {dsub(1:1:numSub).name}';
            
            %display models and subjects
            numModels = numel(models);
            fprintf('%s \n',['Model is: ',models{j}])
            fprintf('%s \n', 'Subjects are:')
            display(SubFolds.name)
            
            %collect data for each subject
            for k = 1 : numSub
                
                %get and store subject
                subject =SubFolds.name{k};
                
                %set this subject path and go there
                thisModelsandSubPath = [thisModelPath,'/',subject];
                cd(thisModelsandSubPath)
                
                %load and list data
                datalisting = dir('*datafit*.mat');
                load(datalisting.name);
                
                %mean Cross-valudated R-squared
                meanCrossValmaxlogl(j,k) = output.CrossVmaxlogl;
            end
        end
        
        
        %position data
        for i = 1 : numModels
            x1(i,:) = [1:numModels+1:(numModels+1)*numSub-numModels]+(i-1);
        end
        
        %colors
        colors=[0 0 0;
            .5 .5 .5;
            .7 .7 .7;
            .8 .8 .8;
            .9 .9 .9];
        
        %status
        fprintf('%s \n','Plotting Cross-validated maxlogl...')
        
        %draw cross validated R-squared
        figure('color','w')
        maxCVmaxlogl = nan(numSub,1);
        for i =  1 : numSub
            b = nan(numModels,1);
            for j = 1 : numModels
                
                %ax
                hold all
                
                %when there is data
                if ~isnan(meanCrossValmaxlogl(j,i))
                    
                    %bars (Stocker and Simoncelli, 2006,NN)
                    b(j) = bar(x1(j,i),meanCrossValmaxlogl(j,i),...
                        0.95,...
                        'edgecolor','none',...
                        'facecolor',colors(j,:),...
                        'displayname',models{j}(7:end));
                else
                    %if missing data
                    scatter(j,0,100,...
                        'k');
                end
                
                %graphics
                maxCVmaxlogl(i) = max(meanCrossValmaxlogl(:));
            end
            
            %graphics
            %labels
            ylabel('Mean cross-validated R-squared','fontsize',14)
            xlabel('Subject','fontsize',14)
        end
        legend(b);
        
        %label subjects
        set(gca,'ytick',(0:0.1:max(maxCVmaxlogl)),'yticklabel',...
            (0:0.1:max(maxCVmaxlogl)))
        set(gca,'xtick',x1(1,:),'xticklabel',mysubjects)
        xlim([0 max(x1(:))+1])
    end
end

%case combined von Mises or bimodal prior
%----------------------------------------
if sum(strcmp(myExperiment,'vonMisesPrior'))==1 && ...
        sum(strcmp(myExperiment,'bimodalPrior'))==1
   
    %---
    %AIC
    %---
    if sum(strcmp(varargs,'AIC')) || sum(strcmp(varargs,'AICcmaes'))
        
        %initialize
        ThisExpForAIC=[];
        TheModels = [];
        TheSub    = [];
        numfitP   = [];
        logL_pertrial = [];
        ThisExperiment =[];
        cohAll  = [];
        pstdAll = [];
        priorModesAll = [];
        dispAll = [];
        
        %Warning
        fprintf('%s \n',['WARNING ! we assumed that von Mises exp. was fitted first...',...
            ' Make sure you fitted vonMises exp. before bimodal exp.'])
        
        %models
        for j = 1 : Folds.nb
            
            %set this model folder
            thisModelPath = strcat([pathFold,'/',Folds.name{j}]);
            
            %get the names and informations about subjects folders found in
            %directory
            %subjects names and info
            for i = 1 : length(subjectsAIC)
                dsub(i) = dir([thisModelPath '/*' num2str(subjectsAIC{i}) '*']);
            end
            
            %get the number of subjects
            numSub = length(dsub);
            
            %get the names of subjects' subfolders
            SubFolds.name = {dsub(1:1:numSub).name}';
            
            %display models and subjects
            fprintf('%s \n', 'Subjects are:')
            fprintf('%s \n', 'Model is:',models{j})
            display(SubFolds.name)
                        
            %subjects' data
            for k = 1 : numSub
                
                %subjects
                subject = SubFolds.name{k};
                
                %subject path
                thisModelsandSubPath = [thisModelPath,'/',subject];
                cd(thisModelsandSubPath)
                
                %load data
                datalisting = dir('*datafit*.mat');
                load(datalisting.name);
                
                %data info 
                numTrials = length(output.Logl_pertrialBestfit);
                TheModels = [TheModels ; repmat({models{j}},numTrials,1)];
                TheSub    = [TheSub ; repmat({subject},numTrials,1)];
                numfitP   = [numfitP ; repmat({sum(~isnan(output.fitP))},numTrials,1)];
                logL_pertrial = [logL_pertrial ; Logl_pertrialBestfit];
                cohAll = [cohAll ; output.StimStrengthVM ; output.StimStrengthBim];
                pstdAll = [pstdAll ; output.pstdVM ; output.pstdBim];
                dispAll = [dispAll ; output.dispVM ; output.dispBim];
                
                %modes
                priorModesBimCell = [];
                for ijk = 1 : size(output.priorModesBim,1)
                    priorModesBimCell{ijk} = output.priorModesBim(ijk,:);
                end
                priorModesVmCell = num2cell(output.priorModesVM);
                priorModesAll = [priorModesAll ; priorModesVmCell; ...
                    priorModesBimCell'];

                %initialize and tag von Mises and bimodal Prior trials
                ThisExperiment = [repmat({'vonMisesTrial'},...
                    numel(output.dataVM),1);...
                    repmat({'bimodalTrial'},numel(output.dataBim),1)];
                ThisExpForAIC = [ThisExpForAIC; ThisExperiment];
            end
        end
        cd(pathFold)
        
        %AICs for conditions of interest (that maximize model comparison)
        %----------------------------------------------------------------
        %status
        fprintf('%s \n','Sorting condition that maximize model comparison...')
        
        %conditions of max difference between competition and Bayes models
        %-----------------------------------------------------------------        
        myCondition = nan(size(logL_pertrial));
        for lmn = 1 : size(logL_pertrial)
        
            %case von mises trials
            %---------------------
            if strcmp(ThisExpForAIC{lmn},'vonMisesTrial')
                
                %distance direction to prior
                Dist2Prior = SLvectors2signedAngle(dispAll(lmn),225,'polar');
        
                %myCondition = abs(Dist2Prior)>=0;  %all data
                myCondition(lmn) = abs(Dist2Prior)>=90; %or distance of at least 90 deg
                
            elseif strcmp(ThisExpForAIC{lmn},'bimodalTrial')
                
                %case bimodal trials
                %-------------------
                %calculate distance of displayed direction to the prior
                %counterclock mode at each trial
                %when the distance between the first and second mode is <0 the
                %first mode is more counterclockwise than the second and when it is
                %>0, the second mode is more clockwise than the first. Now get the
                %coordinates of the counterclock mode (left)
                mode1 = priorModesAll{lmn}(1);
                mode2 = priorModesAll{lmn}(2);
                distModes1Prior = SLvectors2signedAngle(225,mode1,'polar');
                
                %clockwise and counterclockwise modes
                if distModes1Prior > 0
                    cwMode = mode1;
                    ccwMode = mode2;
                else
                    cwMode = mode2;
                    ccwMode = mode1;
                end
                disMeanTocwMode = SLvectors2signedAngle(225,cwMode,'polar');
                disMeanToccwMode = SLvectors2signedAngle(225,ccwMode,'polar');

                %select direction that are not in between the two modes where
                %(3 modes for competition vs. modes for Bayesian). 
                %(tails are 30 deg away from modes)
                disMeantoDisp = SLvectors2signedAngle(225,dispAll(lmn),'polar');
                myCondition(lmn) = disMeantoDisp > disMeanTocwMode || disMeantoDisp < disMeanToccwMode;
            end
        end
          
        %sort trials
        TheModels = TheModels(myCondition==1);
        TheSub    = TheSub(myCondition==1);
        numfitP   = [numfitP{myCondition==1}];
        logL_pertrial = logL_pertrial(myCondition==1);
        cohAll = cohAll(myCondition==1);
        pstdAll = pstdAll(myCondition==1);
        dispAll = dispAll(myCondition==1);
        priorModesAll = priorModesAll(myCondition==1);
        ThisExpForAIC = ThisExpForAIC(myCondition==1);
        
        %AICs
        %----
        %AIC is calculated like this "AIC = 2k - 2ln(L)". We calculate a
        %matrix of AIC values for each subject (rows) and each models
        %(column).
        %ln(L) is the log likelihood of the entire dataset given each
        %model. It is the sum of the log likelihood of each trial-data.
        %ln(L) = ln(L(data1/model1)) + ln(L(data2/model1)) + ....
        %models
        thisModel = unique(TheModels,'stable');
        numModels = length(thisModel);
        
        %last subject
        uniqSubj = unique(TheSub);
        numTotalSub = numel(uniqSubj);
        
        %this model's trials
        %get number of fit parameter
        %get subjects tested for this model
        %SubThisModel=nan(20,numModels);
        %initialize
        AIC = nan(numModels,numTotalSub);
        loglall = nan(size(AIC));
        numTrial = nan(size(AIC));
        numfitPthisModel = nan(numModels,1);
        
        %models
        for i = 1 : numModels
            ThisModelTrials = ismember(TheModels,thisModel{i});
            numfitPthisModel(i) = unique(numfitP(ThisModelTrials),'stable');
            SubThisModel = unique(TheSub(ThisModelTrials));
            
            %subjects
            for j = 1 : numTotalSub
                if sum(strcmp(SubThisModel,uniqSubj{j}))==1
                    
                    %this subject's trials
                    ThisSubTrials = ismember(TheSub,uniqSubj{j});
                    
                    %AIC and DeltaI
                    loglall(i,j) = sum(logL_pertrial(ThisModelTrials & ...
                        ThisSubTrials));
                    AIC(i,j) = 2*numfitPthisModel(i) - 2*loglall(i,j);
                    
                    %number of trials eacn condition
                    numTrial(i,j) = sum(ThisModelTrials & ThisSubTrials);
                else
                    %missing data
                    loglall(i,j)=NaN;
                    AIC(i,j)=NaN;
                end
            end
        end
        
        %draw AICs (model,subjects)
        figure('color','w')
        colors=[0 0 0;
            .5 .5 .5;
            .7 .7 .7;
            .8 .8 .8;
            .9 .9 .9];
        
        %position data
        for i=1:numModels
            x1(i,:)=[1:numModels+1:(numModels+1)*numTotalSub-numModels]+(i-1);
        end
        
        %models and subjects names
        modelnm = unique(TheModels,'stable');
        
        for i=1:numTotalSub
            for j=1:numModels
                hold all
                
                %draw AIC
                %bars (Stocker and Simoncelli, 2006,NN)
                b(j) = bar(x1(j,i),AIC(j,i),...
                    'edgecolor','none',...
                    'facecolor',colors(j,:),...
                    'displayname',modelnm{j});
            end
        end
        xlim([0 x1(end)+1])
        ylabel('AIC','fontsize',14)
        drawPublishAxis
        
        %draw DeltaI (AIC - min AIC)
        figure('color','w')
        hold all
        maxDeltaI = nan(numTotalSub,1);
        SubID = nan(numTotalSub,1);
        
        %subjects
        for i = 1 : numTotalSub
            
            %suject ID
            SubID(i) = str2double(uniqSubj{i}(4:end));
            
            for j = 1 : numModels
                hold all
                
                %when data
                if ~isnan(AIC(j,i))
                    
                    deltaI = AIC(j,i) - min(AIC(:,i));
                    
                    %bars (Stocker and Simoncelli, 2006,NN)
                    %add 1 just for visbility.
                    %same axes
                    space = 0.03;
                    
                    %plot
                    %trick because delta = 0 (best model) is -inf
                    %if logged. Add space to visualize delta = 0.
                    deltaI(deltaI==0) = 1;
                    b(j) = bar(x1(j,i),log10(deltaI) + space,...
                        0.95,...
                        'edgecolor','none',...
                        'facecolor', colors(j,:),...
                        'displayname', modelnm{j}(7:end));
                else
                    %if missing data
                    scatter(j,0,100,'k');
                end
                
                %graphics
                maxDeltaI(i) = max(AIC(:,i) - min(AIC(:,i)));
            end
            
            %labels
            ylabel('Log $\delta$ I','fontsize',14,'interpreter','latex')
            xlabel('Subject','fontsize',14)
        end
        legend(b)
        
        %label subjects
        set(gca,'xtick',x1(1,:),'xticklabel',SubID)
        xlim([0 max(x1(:))+1])
        
        %indicate AIC difference of 2. (significance evidence for the
        %best model)
        plot(linspace(1,max(x1(:))),linspace(log10(2)+space, log10(2)+space),...
            '--k')
        %SLConventionUp(gcf)
    end
end

%clear varargs
clear varargs






%