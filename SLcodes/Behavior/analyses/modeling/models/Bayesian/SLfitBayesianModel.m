

%SLfitBayesianModel.m
%
%
%     author: steeve laquitaine
%       date: 131205 updated 160113
%     status: Complete
%    purpose: Model motion direction estimation data as Bayesian
%             Also fit alternative (suboptimal Bayesian models)
%
%
%      usage:
%
%
%Inputs:
%
%   subjects : e.g., {'sub01'} or {'sub01','sub02'}. Note that when many
%              subjects are inputs data are averaged over subjects not
%              trials.
%      initp : at least 11 models parameters e.g.,
%              [KstrgLLH KmidLLH KweakLLH KveryWeakPrior KweakPrior KStgPrior KveryStgPrior kcardinal Prand Kmotor]
%
%
%
%Other inputs:
%
%
%              'experiment','vonMisesPrior': set von Mises prior data
%                            'bimodalPrior': or bimodal prior data
%                'vonMisesAndBimodalPriors': combined data
%
%                               'MAPReadout': Bayes with MAP readout
%                          'SamplingReadout': Bayes with sampling readout
%
%                             'withCardinal': cardinal prior
%                          'withoutCardinal': no cardinal prior
%
%                    'filename','myFilename': set saved file name
%
%
%Specific analyses:
%
% Bayesian models
% --------------
%                            'FatTailPrior': Bayes with priors with different strengths & same long tails
%              'FatTailPrior','ChangeWTail': Bayes with priors with different strengths & long tails
% 'FatTailPriorAndLLH','FitEachKvmAndTails': Bayes with priors and LLH with different strength and long tails
%                                            17 initp: [KstrgLLH KmidLLH KweakLLH KveryWeakPrior KweakPrior KStgPrior KveryStgPrior kcardinal Prand Kmotor 3Liketails 4PriorTails]
%
%
% fit
% ---
%            'MaxLikelihoodFit','fminsearch': fit models to trial data by minimizing neglogl with fminsearch (max likelihood)
%                 'MaxLikelihoodFit','cmaes': fit models to trial data by minimizing neglogl with cmaes (max likelihood)
%                            'CrossValR2Fit': fit models to data mean and std per conditions (cross-validated R2)
%                                             least-square fit for cross-Validation
%                                             this requires very good estimates of estimates
%                                             mean and std., i.e., lots of sample per condition
%                   'CrossValR2FitTrialData': fit trial data
%                                             R2 is cross-validated
%
%
% Best predictions
% ---------------
%
%                                     initP : a vector of "NaN" to get existing best params
%                                             in 'directoryFitParameter'
% 'directoryFitParameter','../modelfit/AIC/': path of model best fit parameters
%                                             a .mat file with the tag
%                                             "datafit" must already exist
%                                             in the indicated directory.
%
%                                             The file must contain the fit
%                                             parameters in variable "fitP.p"
%
%              'modelPredictions','withData': show predictions for best fit
%                                             parameters  retrieved
%                                             from "directoryFitParameter"
%                                             path (if initP are NaN) or input as
%                                             initp.
%                                   'noData': draw predictions for input
%                                             parameters without showing data
%
%
%    simulation
%    ----------
%
%                                           initP : a vector of of parameters
%
%                   'modelPredictions','withData',: simulated predictions
%                                        ,'noData': draw predictions for input
%                                                   parameters without showing data
%
%                             'inputFitParameters': tag for simulation
%
%
%           Simulate fake data
%           ------------------
%
%                              'testFakedata': to simulate data
%                                              The following variable must
%                                              exist in the workspace :
%                                                      data: series of 1 by Nestimates
%                                                      disp: series of 1 by Nmotion directions
%                                              StimStrength: series of 1 by Ncoherences
%                                                      pstd: series of 1 by Nprior std
%                                                priorModes: series of 1 by Nmodes
%
%
%                               'LoadDataBase': to quickly load recent database "datbank.mat"
%          'dataPath','/Users/steeve/.../proj': to input data path
%
%
%options
%
%                                        debug: maxiter = 1
%
% References:
%     -Hurliman et al, 2002,VR
%     -Stocker&Simoncelli,2006,NN
%     -Girshick&Simoncelli,2011,NN
%     -Chalk&Series,2012,JoV
%
%
%If you need examples do:
%
%     help SLfitBayesianModelExamples.m


function [fitP,fitPbkp,R2,sdata,fitPt,negLogl,negLoglbkp,...
    Logl_pertrialBestfit,output] = SLfitBayesianModel(subjects,...
    initp,...
    varargin)

%time
t0 = tic;

%help
if ieNotDefined('subjects')
    help SLfitBayesianModel
    return
end

%initialize outputs
output = [];

%backup .m file
%Cluster
%doesn't work on stanford sherlock cluster
if ~slIsInput(varargin,'sherlock')
    %     mfilename = SLgetActivemFile;
    %     myMfile = SLbackupMfileInWSpace(mfilename);
    %     output.mfilename  = mfilename;
    %     output.myMfile = myMfile;
end

%data path
%---------
%von Mises prior experiment
vrg = varargin;
%local computer
if ~slIsInput(vrg,'sherlock')
    if slIsInput(vrg,'vonMisesPrior')&&~slIsInput(vrg,'bimodalPrior')
        %get data
        if slIsInput(vrg,'dataPathVM')
            dataPathVM = varargin{find(strcmp(varargin,'dataPathVM'))+1};
            cd(dataPathVM)
        else
            fprintf('%s \n','(SLfitBayesianModel) You need to set dataPath ...')
            dataPathVM ='~/data/dataPsychophy/proj01_priorStrength';%uigetdir(cd,'Pickup your project e.g., /dataPsychophy/Exp01...');
            cd(dataPathVM)
        end
        vararginVM = [varargin,'dataPath',dataPathVM];
        %Load or create database
        if sum(strcmp(varargin,'LoadDataBase'))==0
            fprintf('%s \n','(SLfitBayesianModel) Creating database ...')
            databankVM = SLMakedatabank(subjects,vararginVM);
            databankBim = [];
        else
            %if database exists
            cd data
            if ~isempty(dir('*datbank.mat'))
                fprintf('%s \n','(SLfitBayesianModel) Loading databank in current directory...')
                load('datbank')
                databankVM = databank;
            else
                fprintf('%s \n','(SLfitBayesianModel) Database file datbank.mat does not exist...')
                keyboard
            end
        end
    end
    %cluster
elseif slIsInput(varargin,'sherlock')
    vararginVM  = [varargin,'dataPath','~/data'];
    databankVM  = SLMakedatabank(subjects,vararginVM);
    databankBim = [];
end

%bimodal prior experiment
if slIsInput(varargin,'bimodalPrior') && ~slIsInput(varargin,'vonMisesPrior')
    %get data
    if sum(strcmp(varargin,'dataPathBim'))
        dataPathBim = varargin{find(strcmp(varargin,'dataPathBim'))+1};
        cd(dataPathBim)
    else
        fprintf('%s \n','(SLfitBayesianModel) You need to set dataPath')
        return
    end
    databankVM = [];
    vararginBim = [varargin,'dataPath',dataPathBim];
    databankBim = SLMakedatabank(subjects,vararginBim);
end

%(case combined von Mises and bimodal priors)
%we create two databanks: one for von Mises and one for bimodal priors
if sum(strcmp(varargin,'vonMisesPrior')) && ...
        sum(strcmp(varargin,'bimodalPrior'))
    
    %set,go and get databank von Mises prior data path
    if sum(strcmp(varargin,'dataPathVM'))
        dataPathVM = varargin{find(strcmp(varargin,'dataPathVM'))+1};
        cd(dataPathVM)
    else
        fprintf('%s \n','(SLfitBayesianModel) You need to set dataPathVM')
        return
    end
    vararginVM = [varargin,'dataPath',dataPathVM];
    vararginVM(strcmp(vararginVM,'bimodalPrior'))=[];
    databankVM = SLMakedatabank(subjects,vararginVM);
    
    %set,go and get databank bimodal prior data path.
    if sum(strcmp(varargin,'dataPathBim'))
        dataPathBim = varargin{find(strcmp(varargin,'dataPathBim'))+1};
        cd(dataPathBim)
    else
        fprintf('%s \n','(SLfitBayesianModel) You need to set dataPathBim')
        return
    end
    vararginBim = [varargin,'dataPath',dataPathBim];
    vararginBim(strcmp(vararginBim,'vonMisesPrior'))=[];
    databankBim = SLMakedatabank(subjects,vararginBim);
end


%for prior modes
%(case von Mises prior)
if sum(strcmp(varargin,'experiment'))
    if sum(strcmp(varargin,'vonMisesPrior')) && ...
            sum(strcmp(varargin,'bimodalPrior'))==0
        priorShape = 'vonMisesPrior';
    end
end

%(case bimodal prior)
if sum(strcmp(varargin,'experiment'))
    if sum(strcmp(varargin,'bimodalPrior')) && ...
            sum(strcmp(varargin,'vonMisesPrior'))==0
        priorShape = 'bimodalPrior';
    end
end

%(case combined von Mises and bimodal prior)
if sum(strcmp(varargin,'experiment'))
    if sum(strcmp(varargin,'vonMisesPrior')) && ...
            sum(strcmp(varargin,'bimodalPrior'))
        priorShape = {'vonMisesPrior','bimodalPrior'};
    end
end

%backup file name
if sum(strcmp(varargin,'filename'))
    filename = varargin{find(strcmp(varargin,'filename'))+1};
else
    filename  = input('SLfitBayesianModel) Please set the name of the .mat file that will be saved. (e.g., dataSub01)','s');
end

%check initp
if size(initp) < 7
    error('(SLfitBayesianModel) Please check you initp. Something is wrong...')
    keyboard
end

%(case cardinal prior or not)
if isnan(initp(8))
    TheModel='withoutCardinal';
else
    TheModel='withCardinal';
end

%Initialize outputs
%(case maximum likelihood fit)
R2      =[];
udata   =[];
sdata   =[];
dataOut =[];
predOut =[];
FAbp    =[];
FA      =[];
fitPt   =[];
fitPbkp =[];
negLoglbkp=[];
Logl_pertrialBestfit=[];
negLogl=[];

%(case cross-validated R2)
minSSE = [];
SSE_bkp = [];

%You can set fit parameterss (or not,'[]').
fitP = [];

%check that an analysis is input
posTarget = SLfindword(varargin,{'MaxLikelihoodFit',...
    'CrossValR2Fit',...
    'CrossValR2FitTrialData',...
    'modelPredictions'});

if sum(posTarget)==0
    
    %status
    fprintf('%s \n',['(SLfitBayesianModel) You need to input an analysis:',...
        ' Possible analyses are: '],...
        ' - MaxLikelihoodFit, ',...
        ' - CrossValR2Fit, ',...
        ' - CrossValR2FitTrialData, ',...
        ' - modelPredictions')
end

%(case Maximum likelihood fit)
if sum(strcmp(varargin,'MaxLikelihoodFit'))
    
    %status
    fprintf('%s \n',['(MLft5) Now fitting the model with maximum',...
        ' likelihood method....'])
    
    %fit
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
        MLft5(databankVM,databankBim,fitP,initp,priorShape,filename,...
        output,varargin);
    
    %save data
    output.duration = toc(t0);
    mkdir(['fitBayes' subjects{1}])
    cd(['fitBayes' subjects{1}])
    slPrintfStr('SLfitBayesianModel',['Saving data in' pwd])
    save(filename)
end

%(case Cross-Validated R-squared fit)
if sum(strcmp(varargin,'CrossValR2Fit'))
    
    %status
    fprintf('%s \n',['(CrossValR2fit) Fitting the model for',...
        ' cross-validated R2 method....'])
    
    %Coding is quick for now. Ideally I should implement the 3 conditions
    %'vom Mises' only, 'bimodal only' or 'combined' as I did for MLft5
    %function.
    %case von Mises prior only
    if sum(strcmp(varargin,'vonMisesPrior'))
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2fit(databankVM,fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %case bimodal prior only
    if sum(strcmp(varargin,'bimodalPrior'))
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2fit(databankBim,fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %status
    fprintf('%s \n','(CrossValR2fit) Ending Cross-validation.')
    
    %save data
    save(filename)
    fprintf('%s \n','(CrossValR2fit) .mat file is saved.')
end

%(case Cross-Validated R^2 fit with trial data)
if sum(strcmp(varargin,'CrossValR2FitTrialData'))==1
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialData) Now fitting the model with',...
        'cross-validated R2 method based on single trial data method....'])
    
    %case von Mises prior only
    if sum(strcmp(varargin,'vonMisesPrior'))
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2FitTrialData(databankVM,...
            fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %case bimodal prior only
    if sum(strcmp(varargin,'bimodalPrior'))
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2FitTrialData(databankBim,fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %status
    fprintf('%s \n','(CrossValR2fit) Ending Cross-validation.')
    
    %save data
    save(filename)
    fprintf('%s \n','(CrossValR2fit) .mat file is saved.')
end

%Case the model uses priors with different tail weights and
%the same von Mises strength
if sum(strcmp(varargin,'ChangeWTail'))==1
    fprintf('%s \n',['(SLfitBayesianModel) Your Bayesian model priors will have',...
        ' different tails weights and same von Mises strengths'])
end

%Case priors and LLH with different tail weights and von Mises strengths
if sum(strcmp(varargin,'FitEachKvmAndTails'))==1
    fprintf('%s \n',['(SLfitBayesianModel) WARNING: Make sure fit parameters are in this order:',...
        'kl24, kl12, kl6, klearnt1, klearnt2, klearnt3, klearnt4, kcard, Prand, Kmotor',...
        'Tail_llh1, Tail_llh2, Tail_llh3, Tail_prior1, Tail_prior2, Tail_prior3, Tail_prior4'])
end

%(If we want data mean, std and distribution with models'predictions)
if sum(strcmp(varargin,'modelPredictions'))==1
    
    %case von Mises prior
    if sum(strcmp(varargin,'vonMisesPrior')) && sum(strcmp(varargin,'bimodalPrior'))==0
        databank = databankVM;
    end
    
    %case bimodal prior
    if sum(strcmp(varargin,'bimodalPrior')) && sum(strcmp(varargin,'vonMisesPrior'))==0
        databank = databankBim;
    end
    
    %plot data or not
    if sum(strcmp(varargin,'withData'))==0 && sum(strcmp(varargin,'noData'))==0
        
        %please input with or without data
        error('(SLfitBayesianModel) Please input "withData" or "noData" argument to display or not the data...')
        keyboard
        
    end
    
    %with data
    if sum(strcmp(varargin,'withData'))==1
        
        %status
        fprintf('%s \n',['(SLdrawModelsPredictionCentered) Now draw data mean,'...
            'std and distributions with best models predictions'])
        
        %------
        %models
        %------
        %Bayes Sampling
        if any(strcmp(varargin,'SamplingReadout'))==1 && ...
                sum(strcmp(varargin,'FitEachKvmAndTails'))==0
            
            fprintf('(SLdrawModelsPredictionCentered) Your model readout is Sampling.')
            modelFold = 'model_Bayes_Sampling';
            
            %Bayes MAP
        elseif any(strcmp(varargin,'MAPReadout'))==1 && any(strcmp(varargin,'withCardinal'))==0
            
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Your model readout is MAP ')
            modelFold='model_Bayes_MAP';
            
            %Bayes MAP with card
        elseif any(strcmp(varargin,'MAPReadout'))==1 && any(strcmp(varargin,'withCardinal'))==1
            
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Your model is Bayes MAP with cardinal ')
            modelFold='model_Bayes_MAP_withCard';
            
            %Heavy-tailed LLH and Prior Bayes Sampling (free Tails and strengths)
        elseif sum(strcmp(varargin,'FatTailPriorAndLLH'))==1 && ...
                sum(strcmp(varargin,'FitEachKvmAndTails'))==1 && ...
                sum(strcmp(varargin,'SamplingReadout'))==1
            
            fprintf('(SLdrawModelsPredictionCentered) Your model is Heavy-tailed LLH and Prior Bayes Sampling (free Tails and strengths) \n')
            modelFold='model_Bayes_Sampling_FatTailPriorAndLLH_FitEachKvmAndTails';
            
        else
            fprintf('(SLdrawModelsPredictionCentered) You need to input a readout ...\n')
            keyboard
        end
        
        %number of subjects input
        inputSubjects = subjects;
        subjectsNum = unique(databank.subjects);
        numSub = numel(subjectsNum);
        
        %warn that data and predictions will be averaged across subjects
        if numSub > 1
            fprintf(['(SLdrawModelsPredictionCentered) More that 1 subject have been found',...
                '. Data and models predictions will be averaged across subjects /n'])
        end
        
        %all experimental conditions
        %case von Mises prior
        %--------------------
        if strcmp(priorShape,'vonMisesPrior')
            output.uniqCond = SLuniqpair([databank.Pstd databank.stimStrength ...
                databank.stimFeatureDeg]);
            numcond = size(output.uniqCond,1);
            
            %case bimodal prior
            %------------------
        elseif strcmp(priorShape,'bimodalPrior')
            
            %prior conditions
            priorModes = cell2mat(databank.priormodes);
            priorCond = priorModes(:,2) - priorModes(:,1);
            
            %check that all prior conditions are there
            numPriorCond = size(SLuniqpair(priorModes),1);
            
            if numel(unique(priorCond)) == numPriorCond
                
                output.uniqCond = SLuniqpair([priorCond ...
                    databank.stimStrength ...
                    databank.stimFeatureDeg]);
                numcond = size(output.uniqCond,1);
                
            end
            
            %Warning if priorShape is missing
        else
            fprintf('%s \n',['(SLfitBayesianModel : ModelPredictions)',...
                ' You need to input a prior type'])
            keyboard
        end
        
        %preallocate
        output.meanData = nan(numcond,numSub);
        output.stdData  = nan(numcond,numSub);
        output.meanPred = nan(numcond,numSub);
        output.stdPred  = nan(numcond,numSub);
        
        %get subject's best fit parameters
        %data come from modelfit folder. The tree structure of the folders
        %within the  experiment folder must remain the same. The name of
        %the folder must not change
        %!! can probably improve here by directly getting all data from
        %modelfit folder.
        
        %PARAMETERS
        %----------
        %case directory for best fit parameters is input
        %-----------------------------------------------
        if sum(strcmp(varargin,'directoryFitParameter'))
            
            %directory
            inputdir = varargin{find(strcmp(varargin,...
                'directoryFitParameter')) + 1};
            
            dirModelFit = dir([inputdir,modelFold,'/sub*']);
            
            %Check that the directory exists
            if isempty(dirModelFit)==1
                fprintf(['(SLfitBayesianModel) !! WARNING !! Sorry but',inputdir,modelFold,'/sub does not exist.....\n'])
                keyboard
            end
            
            %subjects in directory
            %subjectsInDir = {dirModelFit.name};
            subjectsInDir = subjects;
            numSubInDir = numel(subjectsInDir);
            
            %status
            fprintf('(SLfitBayesianModel) The best fit parameters path was input manually \n')
            
            %case we use input Parameters
            %----------------------------
        elseif sum(strcmp(varargin,'inputFitParameters'))
            
            inputdir = initp;
            numSubInDir = numSub;
            subjectsInDir = inputSubjects;
            
            %case not input
            %--------------
        else
            error('(SLfitBayesianModel) Please input best fit parameters directory after "directoryFitParameter"...e.g., /modelfit/AIC/')
            keyboard
        end
        
        %list of subjects
        subjectsInDirlist = nan(numSubInDir,1);
        for ij = 1 : numel(subjectsInDir)
            subjectsInDirlist(ij) = str2double(subjectsInDir{ij}(4:end));
        end              
        
        %data and predictions for each subject
        for sub = 1 : numSub
            
            %status
            fprintf(['(SLfitBayesianModel) Retrieving data for subject(s) : ',...
                num2str(subjectsNum(sub)),'\n'])
            
            %this trial's subject
            TrialsThisSub = subjectsNum(sub)==databank.subjects;
            
            %data and conditions
            data = round(databank.estimatesDeg(TrialsThisSub));
            d    = databank.stimFeatureDeg(TrialsThisSub);
            coh  = databank.stimStrength(TrialsThisSub);
            pstd = databank.Pstd(TrialsThisSub);
            priorModes = cell2mat(databank.priormodes(TrialsThisSub,:));
            
            %data mean and std            
            [meanDatatmp,stdDatatmp,dataCondtmp] = SLmakeDataMeanAndStd(data,...
                d,coh,pstd,priorModes,priorShape);
            

            
            
            
            
            
            
            %if we want best fit parameters, show data and predictions
            %---------------------------------------------------------
            %Directory for best fit parameters from maximum likelihood fit
            %must be input (e.g.,.../modelfit/..AIC folder).
            if sum(isnan(initp))==length(initp)
                
                %find subjectsInDir that matches subjects append this subject to
                %loading directory, find datafit mat file, append to directory and
                %load model's best fit parameters
                if sum(strcmp(varargin,'directoryFitParameter'))
                    
                    %subject
                    thisSubInDir = subjectsInDir(subjectsInDirlist==subjectsNum(sub));
                    
                    %file
                    FitfilenameThisSub = dir([inputdir,modelFold,'/' ...
                        cell2mat(thisSubInDir) '/datafit*']);
                    FitfilenameThisSub = FitfilenameThisSub.name;
                    
                    %load
                    load([inputdir,modelFold,'/',thisSubInDir{1},'/',...
                        FitfilenameThisSub(1:end-4)],'fitP')
                    
                    %Warning
                else
                    fprintf('%s \n',['(SLfitBayesianModel) Please input',...
                        'directory for best fit parameters after varargin',...
                        ' "directoryFitParameter"'])
                end
                %                 %case von Mises prior
                %                 %--------------------
                %                 if strcmp(priorShape,'vonMisesPrior')
                %                     thisSubInDir = subjectsInDir(subjectsInDirlist==subjects(sub));
                %                     FitfilenameThisSub = dir(['../modelfit/AIC/',modelFold,'/' ...
                %                         cell2mat(thisSubInDir) '/datafit*']);
                %                     FitfilenameThisSub = FitfilenameThisSub.name;
                %                     load(['../modelfit/AIC/',modelFold,'/' thisSubInDir{1},'/',...
                %                         FitfilenameThisSub(1:end-4)],'fitP')
                %                 end
                %
                %                 %case bimodal prior
                %                 %--------------------
                %                 if strcmp(priorShape,'bimodalPrior')
                %                     thisSubInDir = subjectsInDir(subjectsInDirlist==subjects(sub));
                %                     FitfilenameThisSub = dir(['../modelfit/bimodalPrior/AIC/',modelFold,'/' ...
                %                         cell2mat(thisSubInDir) '/datafit*']);
                %                     FitfilenameThisSub = FitfilenameThisSub.name;
                %                     load(['../modelfit/bimodalPrior/AIC/',modelFold,'/' thisSubInDir{1},'/',...
                %                         FitfilenameThisSub(1:end-4)],'fitP')
                %                 end
                %
                %                 %case combined prior
                %                 %------------------
                %                 if strcmp(priorShape,'combinedPrior')
                %                     thisSubInDir = subjectsInDir(subjectsInDirlist==subjects(sub));
                %                     FitfilenameThisSub = dir(['../modelfit/combinedPrior/AIC/',modelFold,'/' ...
                %                         cell2mat(thisSubInDir) '/datafit*']);
                %                     FitfilenameThisSub = FitfilenameThisSub.name;
                %                     load(['../modelfit/AIC/',modelFold,'/' thisSubInDir{1},'/',...
                %                         FitfilenameThisSub(1:end-4)],'fitP')
                %                 end
                
                %store subjects' best fit parameters
                output.fitP(sub,:) = fitP.p;
                
                %add NaN weightFatTail parameter to old data
                if length(output.fitP(sub,:)) < 11
                    
                    output.fitP(sub,end+1) = NaN;
                    
                    %status
                    fprintf('%s \n',['(SLfitBayesianModel) This is an old set of parameters.',...
                        ' NaN weight Fat Tail prior has been added to',...
                        ' it making a set of 11 para.'])
                    
                end
                
                clear fitP
            end
            
            %case simulation parameters are input (free simulation)
            %------------------------------------------------------
            if sum(isnan(initp))~=length(initp)
                output.fitP(sub,:) = initp;
                
                %status
                fprintf('%s \n','(SLfitBayesianModel) Using input model parameters')
            end
            
            %models' estimate mean, std and distribution predictions
            [meanPredtmp,stdPredtmp,cond,PestimateGivenModelUniq,...
                MAP,output] = SLmakePredictionsBayesianModel(data,d,coh,pstd,output.fitP(sub,:),priorShape,...
                priorModes,TheModel,[],output,varargin);
            
            fprintf('---------------------------------------------------- \n')
            disp(output)
            fprintf('---------------------------------------------------- \n')
            
            %subject's estimate distribution
            [pData,xpdf] = makeDataDist(data,d,coh,pstd,priorModes,cond,...
                priorShape);
            
            %make sure predicted and data distributions are calculated on
            %the same space. Works for predicted estimate space = 1:1:360 deg.
            %adjust predicted estimate probability distribution from 1:1:360 to
            %0:10:360 by summing probabilities within consecutive bins of 10
            %deg (law of probabilities).
            commonSpace = xpdf;
            numSpace = numel(commonSpace)-1;
            [~,bins]= histc(0:1:360,commonSpace);
            bins(end) = [];
            bins = bins';
            numCond = size(cond,1);
            pPred = nan(numSpace,numCond);
            for ijk = 1 : numCond
                pPred(:,ijk) = SLcumSumInBin(PestimateGivenModelUniq(:,ijk),bins);
            end
            commonSpace = commonSpace(2:end);
            
            %match conditions, data and predictions across predictions and
            %data and across subjects
            numPredEst = size(pPred,1);
            numDataEst = size(pData,1);
            for j = 1 : size(output.uniqCond,1)
                
                %case cond (predictions) and dataCondtmp (data) are same as
                %expected
                %get current condition sync with a reference ordering of the
                %conditions
                if isequal(dataCondtmp,cond)
                    [~,posThisCond] = ismember(output.uniqCond(j,:),...
                        dataCondtmp,'rows');
                    
                    %case condition exists
                    if ~isequal(posThisCond,0)
                        
                        %match data and predictions mean, std and distributions
                        %between conditions and subjects to the order in
                        %output.uniqCond
                        output.meanData(j,sub) = meanDatatmp(posThisCond);
                        output.stdData(j,sub)  = stdDatatmp(posThisCond);
                        output.meanPred(j,sub) = meanPredtmp(posThisCond);
                        output.stdPred(j,sub)  = stdPredtmp(posThisCond);
                        output.PdisPred(:,j,sub) = pPred(:,posThisCond);
                        output.PdisData(:,j,sub) = pData(:,posThisCond);
                    else
                        %else NaN
                        output.meanData(j,sub) = NaN;
                        output.stdData(j,sub)  = NaN;
                        output.meanPred(j,sub) = NaN;
                        output.stdPred(j,sub)  = NaN;
                        output.PdisPred(:,j,sub) = NaN(numPredEst,1);
                        output.PdisData(:,j,sub) = NaN(numDataEst,1);
                    end
                end
            end
        end
        
        %stop or continue ?
        YesNo = input('Do you want continue: y/n ?','s');
        if strcmp(YesNo,'n'); keyboard; end
        
                
        
        %---------------------
        %average over subjects
        %---------------------
        
        %calculate data mean, std and distributions
        for i = 1 : size(output.meanData,1)
            CircMeanOvSub = SLcircMeanStd(output.meanData(i,:)','polar');
            output.meanDataOvSub(i) = CircMeanOvSub.deg.mean;
            output.semOfMeanDataOvSub(i) = CircMeanOvSub.deg.sem;
        end
        output.meanDataOvSub    = SLmakeColumn(output.meanDataOvSub);
        output.stdDataOvSub     = nanmean(output.stdData,2);
        output.semOfstdDataOvSub = sem(output.stdData,1)';
        output.meanDisDataOvSub = nanmean(output.PdisData,3);       
        
        %calculate predictions mean and std and distributions
        for i = 1 : size(output.meanPred,1)
            PredCircMeanOvSub = SLcircMeanStd(output.meanPred(i,:)','polar');
            output.meanPredOvSub(i) = PredCircMeanOvSub.deg.mean;
            output.semOfMeanPredOvSub(i) = PredCircMeanOvSub.deg.sem;
        end
        output.meanPredOvSub    = SLmakeColumn(output.meanPredOvSub);
        output.stdPredOvSub     = nanmean(output.stdPred,2);
        output.semOfstdPredOvSub = sem(output.stdPred,1)';
        output.meanDisPredOvSub = nanmean(output.PdisPred,3);
        
        %draw data mean and std and their models' predictions
        SLdrawModelsPredictionCentered(output.meanDataOvSub,...
            output.semOfMeanDataOvSub,output.stdDataOvSub,...
            output.semOfstdDataOvSub,output.meanPredOvSub,...
            output.semOfMeanPredOvSub,output.stdPredOvSub,...
            output.semOfstdPredOvSub,output.uniqCond,priorModes,...
            priorShape,'yCentered')
        SLConventionUp
        
        %draw predictions about estimates' distributions
        SLdrawModelsPredictionHistwithDataCentered(output.meanDisDataOvSub,commonSpace,...
            output.meanDisPredOvSub,d,coh,pstd,priorModes,output.uniqCond,...
            priorShape,varargin)
        
    end
    
    %-----------------------------------
    %Only predictions - do not draw data
    %-----------------------------------
    if sum(strcmp(varargin,'noData'))==1
        
        %Case input param & predictions without data
        if sum(isnan(initp))~=length(initp)
            
            %status
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Simulating model predictions...')
            
            %status: model readout
            if sum(strcmp(varargin,'SamplingReadout'))==1
                fprintf('%s \n','(SLdrawModelsPredictionCentered) Model readout is Sampling ')
            end
            if sum(strcmp(varargin,'MAPReadout'))==1
                fprintf('%s \n','(SLdrawModelsPredictionCentered) Model readout is MAP ')
            end
            if sum(strcmp(varargin,'SamplingReadout'))==0 && sum(strcmp(varargin,'MAPReadout'))==0
                fprintf('%s \n','(SLdrawModelsPredictionCentered) You need to input a readout ...')
                keyboard
            end
            
            %experimental conditions
            d    = databank.stimFeatureDeg;
            pstd = databank.Pstd;
            coh  = databank.stimStrength;
            priorModes = cell2mat(databank.priormodes);
            output.uniqCond = SLuniqpair([pstd coh d]);
            numcond = size(output.uniqCond,1);
            
            %models' estimate mean, std and distribution predictions
            output.fitP = initp;
            [meanPred,stdPred,cond,PestimateGivenModelUniq,...
                MAP] = SLmakePredictionsBayesianModel([],d,...
                coh,...
                pstd,...
                output.fitP,...
                priorShape,...
                priorModes,...
                TheModel,[],...
                output,...
                varargin);
            
            %make sure predicted and data distributions are calculated on
            %the same space. Works for predicted estimate space = 1:1:360 deg.
            %adjust predicted estimate probability distribution from 1:1:360 to
            %0:10:360 by summing probabilities within consecutive bins of 10
            %deg (law of probabilities).
            %commonSpace = 0:10:360;
            commonSpace = 0:10:360;
            numSpace = numel(commonSpace)-1;
            [~,bins] = histc(0:1:360,commonSpace);
            bins(end) = [];
            bins = bins';
            numCond = size(cond,1);
            pPred = nan(numSpace,numCond);
            for ijk = 1 : numCond
                pPred(:,ijk) = SLcumSumInBin(PestimateGivenModelUniq(:,ijk),bins);
            end
            commonSpace = commonSpace(2:end);
            
            %predicted mean,std and distributions
            output.meanPred    = meanPred;
            output.stdPred     = stdPred;
            output.meanDisPred = pPred;
            output.uniqCond    = cond;
            
            %draw predicted mean, std and distribution
            %mean and std
            SLdrawModelsPredictionCentered([],[],output.meanPred,output.stdPred,...
                output.uniqCond,priorModes,priorShape,'yCentered')
            SLConventionUp
            
            %distribution
            SLdrawModelsPredictionHistwithDataCentered([],commonSpace,...
                output.meanDisPred,...
                d,coh,pstd,priorModes,output.uniqCond,priorShape,varargin)
        end
    end
end

%-----------------------
%Case get parameters std
%-----------------------
if sum(strcmp(varargin,'stdBestfitP'))==1
    %status
    fprintf('%s \n',['(MLft5) Now fitting the model with',...
        'maximum likelihood method to get best fit parameters standard error'])
    %fit
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
        MLft5(databankVM,databankBim,fitP,initp,priorShape,filename,...
        output,varargin);
    %save data
    save(filename)
end

%time
output.fitDuration = toc(t0);

%Make sure everything's saved in this directory
if length(subjects)==1
    mkdir(['fitSLfitBayes' subjects{1}])
    cd(['fitSLfitBayes' subjects{1}])
    exp = varargin{find(strcmp(varargin,'experiment'))+1};
    save(['datafit', subjects{1}(end-1:end),'_',exp,'_SLfitBayes'],'output')
    fprintf('%s %s %s \n','(SLfitBayesianModel) Saving fit results in', ['"data/fitSLfitBayes',subjects{1}],'"')
end







%------------------
%MAX LIKELIHOOD FIT
%------------------
function [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
    MLft5(databankVM,databankBim,fitPinput,...
    initp,...
    priorShape,...
    filename,...
    output,...
    varargin)

%status
fprintf('%s \n','(MLft5) Fitting trial-data with maximum likelihood method \n')

%(case check fitting fake data).
%Get simulated data present in the workspace.
%-------------------------------------------
if sum(strcmp(varargin{1},'testFakedata'))
    
    fprintf('%s \n','------------------------------------------------------------------------------------------------- \n')
    fprintf('%s \n','(MLft5) I am loading your simulated data to check if max likelihood fitting method works properly \n')
    fprintf('%s \n','------------------------------------------------------------------------------------------------- \n')
    
    data = evalin('base','pred');
    disp = evalin('base','d');
    coh = evalin('base','coh');
    pstd = evalin('base','pstd');
    priorModes = evalin('base','priorModes');
    
end

%databank(s)
%-----------
%(case von Mises or bimodal prior)
%---------------------------------
if isempty(databankVM) || isempty(databankBim)
    
    fprintf('%s \n',['(MLft5) Fitting trial-data from the von Mises or ',...
        ' the bimodal prior experiment ... \n'])
    
    %databank
    if ~isempty(databankVM)
        databank = databankVM;
    elseif isempty(databankBim)
        databank = databankBim;
    end
    
    %subjects' data
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'estimatedFeature'))==1)));
    
    %task conditions
    disp = cell2mat(databank.data(:,(strcmp(databank.nm,'FeatureSample'))==1));
    StimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd = cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make 360 and 0 degrees same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    StimStrength  = StimStrength(isnan(data)==0);
    pstd = pstd(isnan(data)==0);
    priorModes = priorModes(isnan(data)==0,:);
    
    %store
    output.data = data;
    output.disp = disp;
    output.StimStrength = StimStrength;
    output.pstd = pstd;
    output.priorModes = priorModes;
end

%(case combined von Mises and bimodal prior)
%-------------------------------------------
if ~isempty(databankVM) && ~isempty(databankBim)
    
    %status
    fprintf('%s \n',['(MLft5) I am fitting the estimates data from combined von Mises and',...
        ' bimodal prior experiments... \n'])
    
    
    %(case von Mises prior)
    %----------------------
    %status
    fprintf('%s \n','(MLft5) I am gathering the von Mises prior estimates data')
    
    %subjects' data
    dataVM = round(cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'estimatedFeature'))==1)));
    
    %experimental conditions
    dispVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'FeatureSample'))==1));
    StimStrengthVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'StimStrength'))==1));
    pstdVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'Pstd'))==1));
    priorModesVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    dataVM(dataVM==0) = 360;
    
    %remove missing data
    dataVM = dataVM(isnan(dataVM)==0);
    dispVM = dispVM(isnan(dataVM)==0);
    StimStrengthVM  = StimStrengthVM(isnan(dataVM)==0);
    pstdVM = pstdVM(isnan(dataVM)==0);
    priorModesVM = priorModesVM(isnan(dataVM)==0,:);
    
    %store
    output.dataVM = dataVM;
    output.dispVM = dispVM;
    output.StimStrengthVM = StimStrengthVM;
    output.pstdVM = pstdVM;
    output.priorModesVM = priorModesVM;
    
    
    %(case bimodal prior)
    %--------------------
    %status
    fprintf('%s \n','(MLft5) I am gathering the bimodal prior estimates data')
    
    %subjects' data
    dataBim = round(cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'estimatedFeature'))==1)));
    
    %experimental conditions
    dispBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'FeatureSample'))==1));
    StimStrengthBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'StimStrength'))==1));
    pstdBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'Pstd'))==1));
    priorModesBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    dataBim(dataBim==0) = 360;
    
    %remove missing data
    dataBim = dataBim(isnan(dataBim)==0);
    dispBim = dispBim(isnan(dataBim)==0);
    StimStrengthBim  = StimStrengthBim(isnan(dataBim)==0);
    pstdBim = pstdBim(isnan(dataBim)==0);
    priorModesBim = priorModesBim(isnan(dataBim)==0,:);
    
    %store
    output.dataBim = dataBim;
    output.dispBim = dispBim;
    output.StimStrengthBim = StimStrengthBim;
    output.pstdBim = pstdBim;
    output.priorModesBim = priorModesBim;
end

%for fminsearch - Nelder-Mead
% OptionsNM=optimset('Display','iter');
OptionsNM = [];

%fit if no input parameters
if isempty(fitPinput)==1
    
    %Sets of initial parameters: best matching parameters are found by
    %matching raw data with simulations graphically (called "Kreal").
    %We pre-generate 27 reasonable sets of initial parameters.
    %+-----+-------+----------------------+
    %|     |       |          prior       |
    %.-----|-------+-------+-------+------+
    %|     |       | true  | strong| weak |
    %.-----|-------+-------+-------+------+
    %|     |true   |(1)t-t |(2)t-s |(3)t-w|
    %.     |-------+-------+-------+------+
    %| llh |strong |(4)s-t |(5)s-s |(6)s-w|
    %|     |-------+-------+-------+------+
    %|     |weak   |(7)w-t |(8)w-s |(9)w-w|
    %+-----+-----------------------+------+
    
    %set 10: likelihood and prior strengths are all same.
    %eventually we may use later
    %10 sets.
    
    %note: strong priors (and weak priors) are 8 times stronger than best
    %matching priors. 8x is the factor for which I see clear deviation
    %of simulation from data.
    
    %We used one intial value for probability of random estimation, motor
    %noise and cardinal prior strength, that we think are relatively small
    %values (high k for motor noise is low motor noise) by looking at the
    %data.
    
    
    %Input initial fit parameters (typically best matching initial values)
    %check that all initial parameters are input
    
    %(case fat tail priors and llh with free vm strengths and tails)
    if sum(strcmp(varargin{:},'FitEachKvmAndTails'))==1
        if length(initp)~=17
            
            %status
            fprintf('%s \n',['(MLft5) Some model initial parameters are missing...',...
                'parameters for Bayesian model with fat tail priors and llh with free vm strengths and tails are:'],...
                '3 kllh', ...
                '4 klearnt',...
                '1 kcardinal',...
                '1 fractRand',...
                '1 motorNoise',...
                '1 weighTail')
            keyboard
        end
        
        %standard Bayesian model
    else
        if length(initp)~=11
            
            %status
            fprintf('%s \n',['(MLft5) Some initial parameters are missing...',...
                'parameters for standard Bayesian model is:', '\n 3 kllh', '\n 4 klearnt','\n 1 kcardinal',...
                '\n 1 fractRand','\n 1 motorNoise','\n 1 weighTail'])
            keyboard
            
        end
        
        kllh       = initp(1:3);
        klearnt    = initp(4:7);
        kcardinal  = initp(8);
        fractRand  = initp(9);
        motorNoise = initp(10);
        weightTail = initp(11);
    end
    
    
    %(Case max LLH fit, Initial para)
    %----------------------------------------------------------------
    if sum(strcmp(varargin{1},'MaxLikelihoodFit'))==1
        
        %if priors are fixed and uniform (no prior)
        if any(strcmp(varargin{:},'fixedPriors')) && all(klearnt)==0
            
            %(1) Input parameters (best matching parameters)
            k0(1,:) = [kllh    klearnt kcardinal fractRand motorNoise weightTail];
            %(4) 8x stronger likelihoods & best matching priors
            k0(2,:) = [kllh.*8 klearnt kcardinal fractRand motorNoise weightTail];
            %(7) 8x weaker likelihood & true priors
            k0(3,:) = [kllh./8 klearnt kcardinal fractRand motorNoise weightTail];
            
            output.initp = k0;
            
            %if priors are free and non-uniform
        else
            %(1) Input parameters (best matching parameters)
            k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise weightTail];
            %(2) 8x stronger priors & best matching llh
            k0(2,:) = [kllh    8*klearnt    kcardinal fractRand motorNoise weightTail];
            %(3) 8x weaker priors & best matching llh
            k0(3,:) = [kllh    klearnt./8   kcardinal fractRand motorNoise weightTail];
            %(4) 8x stronger likelihoods & best matching priors
            k0(4,:) = [kllh.*8 klearnt      kcardinal fractRand motorNoise weightTail];
            %(5) 8x stronger likelihoods & stronger priors
            k0(5,:) = [[kllh   klearnt].*8  kcardinal fractRand motorNoise weightTail];
            %(6) 8x stronger likelihoods & weaker priors
            k0(6,:) = [kllh.*8 klearnt./8   kcardinal fractRand motorNoise weightTail];
            %(7) 8x weaker likelihood & true priors
            k0(7,:) = [kllh./8 klearnt      kcardinal fractRand motorNoise weightTail];
            %(8) 8x weaker likelihood & stronger priors
            k0(8,:) = [kllh./8 klearnt.*8   kcardinal fractRand motorNoise weightTail];
            %(9) 8x weaker likelihood & weaker priors
            k0(9,:) = [[kllh   klearnt]./8  kcardinal fractRand motorNoise weightTail];
            %(10)likelihoods & priors are same.
            k0(10,:) = [repmat(nanmean(k0(1,1:7)),1,7) kcardinal fractRand motorNoise weightTail];
        end
                
        %case priors have heavy Tails with same tail weight and different von Mises
        %strengths
        if ~any(strcmp(varargin{:},'fixedPriors'))
            if sum(strcmp(varargin{1},'FitEachKvmAndTails'))==0
                if sum(strcmp(varargin{:},'ChangeWTail'))==0
                    
                    %(1) Input parameters (best matching parameters)
                    k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise weightTail];
                    %(2) 8x stronger priors & best matching llh
                    k0(2,:) = [kllh    8*klearnt    kcardinal fractRand motorNoise weightTail];
                    %(3) 8x weaker priors & best matching llh
                    k0(3,:) = [kllh    klearnt./8   kcardinal fractRand motorNoise weightTail];
                    %(4) 8x stronger likelihoods & best matching priors
                    k0(4,:) = [kllh.*8 klearnt      kcardinal fractRand motorNoise weightTail];
                    %(5) 8x stronger likelihoods & stronger priors
                    k0(5,:) = [[kllh   klearnt].*8  kcardinal fractRand motorNoise weightTail];
                    %(6) 8x stronger likelihoods & weaker priors
                    k0(6,:) = [kllh.*8 klearnt./8   kcardinal fractRand motorNoise weightTail];
                    %(7) 8x weaker likelihood & true priors
                    k0(7,:) = [kllh./8 klearnt      kcardinal fractRand motorNoise weightTail];
                    %(8) 8x weaker likelihood & stronger priors
                    k0(8,:) = [kllh./8 klearnt.*8   kcardinal fractRand motorNoise weightTail];
                    %(9) 8x weaker likelihood & weaker priors
                    k0(9,:) = [[kllh   klearnt]./8  kcardinal fractRand motorNoise weightTail];
                    %(10)likelihoods & priors are same.
                    k0(10,:) = [repmat(nanmean(k0(1,1:7)),1,7) kcardinal fractRand motorNoise weightTail];
                    
                end
            end
        end

        %case priors have heavy Tails with different tail weight and same von Mises
        %strengths
        if sum(strcmp(varargin{:},'ChangeWTail'))==1
            
            %note: here klearnt are the four tails weights and weightTail
            %is the four von Mises strengths
            %(1) Input parameters (best matching parameters)
            k0(1,:) = [kllh klearnt kcardinal fractRand motorNoise weightTail];
            %(2) 8x stronger priors & best matching llh
            k0(2,:) = [kllh klearnt kcardinal fractRand motorNoise 8*weightTail];
            %(3) 8x weaker priors & best matching llh
            k0(3,:) = [kllh klearnt kcardinal fractRand motorNoise weightTail./8];
            %(4) 8x stronger likelihoods & best matching priors
            k0(4,:) = [kllh.*8 klearnt  kcardinal fractRand motorNoise weightTail];
            %(5) 8x stronger likelihoods & stronger priors
            k0(5,:) = [kllh.*8  klearnt kcardinal fractRand motorNoise weightTail.*8];
            %(6) 8x stronger likelihoods & weaker priors
            k0(6,:) = [kllh.*8 klearnt kcardinal fractRand motorNoise weightTail./8];
            %(7) 8x weaker likelihood & true priors
            k0(7,:) = [kllh./8 klearnt kcardinal fractRand motorNoise weightTail];
            %(8) 8x weaker likelihood & stronger priors
            k0(8,:) = [kllh./8 klearnt kcardinal fractRand motorNoise weightTail.*8];
            %(9) 8x weaker likelihood & weaker priors
            k0(9,:) = [kllh./8 klearnt kcardinal fractRand motorNoise weightTail./8];
            %(10)likelihoods & priors are same.
            k0(10,:) = [repmat(nanmean(k0(1,[1,2,3,8,11])),1,3) klearnt kcardinal ...
                fractRand motorNoise nanmean(k0(1,[1,2,3,8,11]))];
        end
        
        
        %(Case max LLH fit, 'FatTailPriorAndLLH', 'FitEachKvmAndTails', Initial para)
        %priors and llh have fit para. for each vm strengths and tails
        if sum(strcmp(varargin{1},'FatTailPriorAndLLH'))==1
            if sum(strcmp(varargin{1},'FitEachKvmAndTails'))==1
                
                kllh       = initp(1:3);
                klearnt    = initp(4:7);
                kcardinal  = initp(8);
                fractRand  = initp(9);
                motorNoise = initp(10);
                Tailllh    = initp(11:13);
                Tailpriors = initp(14:17);
                
                %(1) Input parameters (best matching parameters)
                k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise Tailllh Tailpriors];
                %(2) 8x stronger priors & best matching llh
                k0(2,:) = [kllh    8*klearnt    kcardinal fractRand motorNoise Tailllh Tailpriors];
                %(3) 8x weaker priors & best matching llh
                k0(3,:) = [kllh    klearnt./8   kcardinal fractRand motorNoise Tailllh Tailpriors];
                %(4) 8x stronger likelihoods & best matching priors
                k0(4,:) = [kllh.*8 klearnt      kcardinal fractRand motorNoise Tailllh Tailpriors];
                %(5) 8x stronger likelihoods & stronger priors
                k0(5,:) = [[kllh   klearnt].*8  kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(6) 8x stronger likelihoods & weaker priors
                k0(6,:) = [kllh.*8 klearnt./8   kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(7) 8x weaker likelihood & true priors
                k0(7,:) = [kllh./8 klearnt      kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(8) 8x weaker likelihood & stronger priors
                k0(8,:) = [kllh./8 klearnt.*8   kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(9) 8x weaker likelihood & weaker priors
                k0(9,:) = [[kllh   klearnt]./8  kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(10)likelihoods & priors are same.
                k0(10,:) = [repmat(nanmean(k0(1,1:7)),1,7) kcardinal fractRand motorNoise  Tailllh Tailpriors];
                %(11)likelihoods & priors are same.
                k0(11,:) = [repmat(nanmean(k0(1,1:7)),1,7) kcardinal fractRand motorNoise  repmat(nanmean([Tailllh Tailpriors]),1,7)];
                
            end
        end
    end
    
    
    %(Case we want to get best fit parameters' standard error, fmincon)
    %------------------------------------------------------------------
    %Only keep input parameters (fit parameters e.g., from a previous fit
    %for examples)
    if sum(strcmp(varargin{1},'stdBestfitP'))==1        
        %status
        fprintf('%s \n',['(MLft5) Setting your input initial parameters as'...
            'the only set of initial parameter'])
        k0(1,:) = [kllh klearnt kcardinal fractRand motorNoise weightTail];
    end
    
    %(case we do not want to fit a cardinal prior).
    %Kcardinal=NaN.
    if isnan(kcardinal)
        TheModel = 'withoutCardinal';
    else
        TheModel = 'withCardinal';
    end
    
    %fitting
    %k can be anything >=0;
    fitPbkp = nan(size(k0,1),size(k0,2));
    negLoglbkp = nan(size(k0,1),1);
    
    %(Case we fit data with maximum likelihood
    %-----------------------------------------
    if sum(strcmp(varargin{1},'MaxLikelihoodFit'))==1
        
        %(case von Mises or Bimodal prior)
        %---------------------------------
        if isempty(databankVM) || isempty(databankBim)
            
            %fitting options
            %---------------
            %quick test
            %options = optimset('MaxIter',1,'MaxFunEvals',1);
            
            %actual run
            %options = optimset('MaxIter',200*size(k0,1),...
            %    'MaxFunEvals',10*200*size(k0,1));
            
            %only for 'FatTailPriorAndLLH','FitEachKvmAndTails' model
            %The default TolFun and TolX are very conservative so we keep
            %them (although it means that fitting will take time)
            %options = optimset('MaxIter',10*200*size(k0,1),...
            %    'MaxFunEvals',10*10*200*size(k0,1));
            
            %==========
            %fminsearch
            %==========
            if any(strcmp(varargin{:},'fminsearch'))
                                
                %prior strengths are free
                if ~any(strcmp(varargin{:},'fixedPriors'))
                    [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp] = slfitBayesfminSlogl(k0,data,disp,StimStrength,...
                        pstd,priorShape,priorModes,TheModel,options,varargin{:});
                    output.freePriors = 1;
                    output.fixedPriors = 0;
                else
                    %prior srengths are fixed (10 funEval per iter:
                    %fit already converges (exitflag = 1) tested
                    %with sub01 - 20 min fit
                    options = optimset('MaxIter',1000,'MaxFunEvals',10000);
                    %fit
                    [fitPbkp,negLoglbkp,exitflagbkp,outputFitbkp] = slfitBayesfminSloglPriorsFixed(k0,data,disp,StimStrength,...
                        pstd,priorShape,priorModes,TheModel,options,varargin{:});
                    output.freePriors = 0;
                    output.fixedPriors = 1;                    
                end
                
                %=====
                %CMAES
                %=====
            elseif any(strcmp(varargin{:},'cmaes'))
                
                %options
                options = cmaes;
                options.LBounds = 0;   %no < 0
                options.TolFun = 1e-1; %stop when obj(t2) - obj(t1) < 1e-4
                
                %debug mode
                if slIsInput(varargin{:},'debug')
                    options.MaxIter = 1; %for debugging
                else
                    %testing with defined maxiter for speed----
                    options.MaxIter = 200*10;
                end
                options.MaxFunEvals = options.MaxIter*10;
                
                %have to set those off otherwise parfor does not work
                options.DispFinal = 0;
                options.DispModulo = 0;
                options.SaveVariables = 0;
                options.SaveFilename = 0;
                options.LogModulo = 0;
                options.LogTime = 0;
                options.LogFilenamePrefix = 0;
                options.LogPlot = 0;
                
                %fit
                parfor i = 1:size(k0,1)
                    
                    %status
                    fprintf('\n',['(MaxLikelihoodFit)',' Setting cmaes parameters : ',num2str(i),'/',num2str(size(k0,1))])
                    
                    %fit
                    [fitPtmp,negLogl,~,stopflag,outputFit] = fitCMAE(data,disp,StimStrength,pstd,k0(i,:),priorShape,priorModes,TheModel,options,varargin{1});
                    
                    %Fit parameters and SSE
                    fitPbkp(i,:) = fitPtmp;
                    negLoglbkp(i) = negLogl;
                    exitflagbkp{i} = stopflag;
                    outputFitbkp{i} = outputFit;
                end
                output.fitAlgo = 'cmaes';
                output.fitoptions = options;
            end
            
            %backup
            output.fitPbkp = fitPbkp;
            output.negLoglbkp = negLoglbkp;
            output.outputFitbkp = outputFitbkp;
            output.exitflagbkp  = exitflagbkp;
            
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
                
                %prior strengths are free
                if ~any(strcmp(varargin{:},'fixedPriors'))
                    [~,~,Logl_pertrialBestfit] = SLgetLoglBayesianModel(data,disp,StimStrength,pstd,fitP.p,...
                        priorShape,priorModes,TheModel,varargin{1});
                else
                    [~,~,Logl_pertrialBestfit] = SLgetLoglBayesianModelUnifPrior(data,disp,StimStrength,pstd,fitP.p,...
                        priorShape,priorModes,TheModel,varargin{1});
                end
                
            elseif any(strcmp(varargin{:},'cmaes'))
                %get rid of NaN
                if strcmp(TheModel,'withoutCardinal') || isnan(fitP.p(8))
                    fitP.p(8) = [];  %tail
                    fitP.p(10) = []; %cardinal
                    TheModel = 'withoutCardinal';
                elseif strcmp(TheModel,'withCardinal') || ~isnan(fitP.p(8))
                    fitP.p(11) = [];  %taik
                end
                [~,~,Logl_pertrialBestfit] = SLgetLoglBayesianModelcmae(fitP.p',data,disp,StimStrength,pstd,...
                    priorShape,priorModes,TheModel,varargin{1});
                %rearrange
                if strcmp(TheModel,'withoutCardinal') || isnan(fitP.p(8))
                    fitPdum = nan(11,1);%card
                    fitPdum(1:7) = fitP.p(1:7);
                    fitPdum(9:10) = fitP.p(8:9);
                    fitP.p = fitPdum;
                end
                fitP.p(11) = nan;
            end
            %backup
            output.Logl_pertrialBestfit = Logl_pertrialBestfit;
            output.fitP = fitP.p;
        end
        
        %(case combined von Mises & Bimodal priors)
        %-----------------------------------------
        if ~isempty(databankVM) && ~isempty(databankBim)
            
            %fitting options
            %---------------
            %(debugging)
            options = optimset('MaxIter',1000,'MaxFunEvals',1000);
            %options = [];
            
            %fit
            parfor i = 1 : size(k0,1)
                %for i = 1 : size(k0,1)
                
                %status
                fprintf('%s \n',['(MLft5)',' Set of initial parmeters: ',...
                    num2str(i),'/',num2str(size(k0,1))])
                
                %Nelder-Mead
                [fitPtmp,negLogl,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                    getLoglCombinedPriors(dataVM,dataBim,...
                    dispVM,dispBim,...
                    StimStrengthVM ,StimStrengthBim,...
                    pstdVM,pstdBim,...
                    fitPtmp,...
                    priorModesVM,priorModesBim,...
                    TheModel,varargin{1}),...
                    k0(i,:),...
                    options);
                
                %Fit parameters, SSE, and fit info
                fitPbkp(i,:)    = fitPtmp;
                negLoglbkp(i)   = negLogl;
                exitflagbkp{i}  = exitflag;
                outputFitbkp{i} = outputFit;
            end
            
            %store
            output.fitPbkp = fitPbkp;
            output.negLoglbkp = negLoglbkp;
            output.outputFitbkp = outputFitbkp;
            output.exitflagbkp  = exitflagbkp;
            
            %Max likelihood.
            [negLogl,position] = min(negLoglbkp(:));
            i = ind2sub(size(negLoglbkp),position);
            
            %optimal parameters and Logl_pertrial
            fitP.p = fitPbkp(i,:);
            [~,~,Logl_pertrialBestfit] = getLoglCombinedPriors(dataVM,...
                dataBim,...
                dispVM,dispBim,...
                StimStrengthVM ,StimStrengthBim,...
                pstdVM,pstdBim,...
                fitP.p,...
                priorModesVM,priorModesBim,...
                TheModel,varargin{1});
            
            %store
            output.Logl_pertrialBestfit = Logl_pertrialBestfit;
            output.fitP = fitP.p;
        end
    end
    
    %(Case we ask for best fit parameters' std
    %-----------------------------------------
    if sum(strcmp(varargin{1},'stdBestfitP'))==1
        
        %check cardinal
        if isnan(kcardinal)
            higherConstKcard = 0;
        else
            higherConstKcard = +inf;
        end
        
        %fit
        %options
        %debugging
        %options = optimset('Display','off','Algorithm','active-set',...
        %    'MaxIter',1,'MaxFunEvals',1);
        options = optimset('Display','on','Algorithm','interior-point');
        
        parfor i = 1 : size(k0,1)
            %for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n','(MLft5) Set of the initial parameters...')
            
            %fmincon
            [fitPtmp,negLogl,exitflag,outputFit,LAMBDA,...
                GRAD,fitHessian] = fmincon( @(fitPtmp) ...
                SLgetLoglBayesianModel(data,...
                disp,...
                StimStrength,...
                pstd,...
                fitPtmp,...
                priorShape,...
                priorModes,...
                TheModel,varargin{1}),...
                k0(i,:),...
                [],[],[],[],...
                [0       0    0    0    0    0    0                0 0    0 0],...
                [3000 3000 3000 3000 3000 3000 3000 higherConstKcard 1 3000 1],...
                [],...
                options);
            
            %Fit parameters and negLogl
            fitPbkp(i,:) = fitPtmp;
            negLoglbkp(i) = negLogl;
            exitflagbkp{i} = exitflag;
            outputFitbkp{i} = outputFit;
            LAMBDAbkp{i} = LAMBDA;
            GRADbkp{i} = GRAD;
            fitHessianbkp{i} = fitHessian;
        end
        
        %store
        output.fitPbkp    = fitPbkp;
        output.negLoglbkp = negLoglbkp;
        output.fitHessianbkp = fitHessianbkp;
        output.exitflagbkp = exitflagbkp;
        output.outputFitbkp = outputFitbkp;
        output.LAMBDAbkp    = LAMBDAbkp;
        output.GRADbkp      = GRADbkp;
        
        %Max likelihood
        [negLogl,position] = min(negLoglbkp(:));
        i = ind2sub(size(negLoglbkp),position);
        
        %optimal parameters and Logl_pertrial
        fitP.p = fitPbkp(i,:);
        output.fitHessian = fitHessianbkp{i};
        
        [~,~,Logl_pertrialBestfit] = SLgetLoglBayesianModel(data,disp,StimStrength,pstd,fitP.p,...
            priorShape,priorModes,TheModel,varargin{1});
    end
    
    %(case von Mises prior)
    %---------------------
    if sum(strcmp(priorShape,'vonMisesPrior')) && ...
            ~sum(strcmp(priorShape,'bimodalPrior'))
        
        %case prior have same tail weight, different von Mises strengths
        if sum(strcmp(varargin{:},'ChangeWTail')==0)
            if sum(strcmp(varargin{:},'FitEachKvmAndTails')==0)
                
                fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','pstd80','pstd40','pstd20','pstd10',...
                    'kcardinal','Prand','km','weightTail'};
            end
        end
        %case prior have different tail weight, same von Mises strengths
        if sum(strcmp(varargin{:},'ChangeWTail')==1)
            
            fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','weightTail80','weightTail40','weightTail20','weightTail10',...
                'kcardinal','Prand','km','Kpriors'};
        end
        
        %case prior and llh have free vm strengths and tail weights
        if sum(strcmp(varargin{:},'FitEachKvmAndTails')==1)
            
            fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','Klearnt80','Klearnt40',...
                'Klearnt20','Klearnt10','kcardinal','Prand','km',...
                'Tailllh24','Tailllh12','Tailllh6','Tailprior80',...
                'Tailprior40','Tailprior20','Tailprior10'};
        end
    end
    
    %(case bimodal prior)
    %---------------------
    if sum(strcmp(priorShape,'bimodalPrior')) && ...
            ~sum(strcmp(priorShape,'vonMisesPrior'))
        fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','pstd145_305','pstd165_285',...
            'pstd185_265','pstd205_245',...
            'kcardinal','Prand','km','weightTail'};
    end
    
    %(case combined von Mises and bimodal prior data)
    %------------------------------------------------
    if sum(strcmp(priorShape,'vonMisesPrior')) && ...
            sum(strcmp(priorShape,'bimodalPrior'))
        
        fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','pstd80modes225or20modes145_305',...
            'pstd40modes225or20modes165_285',...
            'pstd20modes225or20modes185_265',...
            'pstd10modes225or20modes205_245',...
            'kcardinal','Prand','km','weightTail'};
    end
    
    %(case no cardinal prior)
    %------------------------
    %Kcardinal=NaN.
    if isnan(kcardinal)
        fitP.p(strcmp(fitP.nm,'kcardinal')) = NaN;
    end
    
    %parameters std
    if sum(strcmp(varargin{1},'stdBestfitP'))==1
        
        [~,~,~,~,~,output] = SLmakePredictionsBayesianModel([],disp,StimStrength,pstd,fitP.p,priorShape,...
            priorModes,TheModel,'Trial',output,varargin{1});
        output.stdPa = SLgetFitParamStd(output.fitHessian,data,output.TrialPred,...
            fitP.p);
        
    elseif sum(strcmp(varargin{1},'MaxLikelihoodFit'))==1
        
        %simple search
        output.TrialPred = [];
        output.stdPa  = [];
    end
    
    %save data (temporary)
    save(filename)
    
    %if model fit parameters are input manually
elseif isempty(fitPinput)==0
    
    %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
    if isnan(fitPinput(8))
        TheModel = 'withoutCardinal';
    else
        TheModel = 'withCardinal';
    end
    
    [negLogl,~,Logl_pertrialBestfit] = SLgetLoglBayesianModel(data,disp,StimStrength,pstd,fitPinput,...
        priorShape,priorModes,TheModel,varargin{1});
    negLoglbkp = negLogl;
    fitPbkp = fitPinput;
    fitP.p = fitPinput;
    fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','klearnt80','klearnt40','klearnt20',...
        'klearnt10','Kcardinal','Prand','km'};
    
    %store data
    save(filename)
end

%logl for combined priors
function [negLogl,fitP,Logl_pertrial] = getLoglCombinedPriors(dataVM,dataBim,...
    displVM,displBim,StimStrengthVM,StimStrengthBim,pstdVM,pstdBim,fitP,...
    priorModesVM,priorModesBim,TheModel,...
    varargin)

%we use the same fit parameters to get log likelihood for von Mises prior
%data and bimodal prior data separately. We then calculated the logl of the
%overall data set by summing the two.

ticLogL = tic;

%negLogl
%von Mises prior
[negLoglVM,fitP,Logl_pertrialVM] = SLgetLoglBayesianModel(dataVM,displVM,StimStrengthVM,pstdVM,...
    fitP,'vonMisesPrior',priorModesVM,TheModel,varargin{1});

if isempty(Logl_pertrialVM)
    displ('empty')
end

%bimodal prior
fitP(4:7) = fitP(6);
[negLoglBim,fitP,Logl_pertrialBim] = SLgetLoglBayesianModel(dataBim,displBim,StimStrengthBim,...
    pstdBim,fitP,'bimodalPrior',priorModesBim,TheModel,varargin{1});

%combine
negLogl = sum([negLoglVM;negLoglBim]);
Logl_pertrial = [Logl_pertrialVM; Logl_pertrialBim];


%Look at fitting. It is 3X faster without drawing.
ti = toc(ticLogL);
fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f   %.2f\n',...
    negLogl,fitP(1),fitP(2),fitP(3),fitP(4),fitP(5),fitP(6),fitP(7),...
    fitP(8),fitP(9),fitP(10),fitP(11),ti)


%-------------------
%CROSS-VALIDATED FIT
%-------------------
%cross-validated R^2 fit (estimates mean and std)
function [pred,data,disp,StimStrength,pstd,stdPa,output] = CrossValR2fit(databank,...
    fitPinput,...
    initp,...
    priorShape,...
    filename,...
    output,...
    varargin)

%STEPS of the procedure
%1. First, sort data as training and test.
%2. Second, get the best fit parameters that minimize sum of square error
%between training data and predictions.
%3. Use the fit Parameters to generate predictions for the test data set
%and calculate test R2.
%4. Start a new training data set and repeat the procedure.
%5. Store the R2 for each of the 5 sets and calculate the mean R2.

%status
fprintf('%s \n',['(CrossValR2fit) Fitting model to data mean and std',...
    ' for each experimental condition'])

%(case check fitting fake data). Get simulated data in the workspace.
if sum(strcmp(varargin{1},'testFakedata'))
    data = evalin('base','pred');
    disp = evalin('base','d');
    StimStrength = evalin('base','StimStrength');
    pstd = evalin('base','pstd');
    priorModes = evalin('base','priorModes');
    
else
    
    %subjects' data and conditions
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'estimatedFeature'))==1)));
    disp = cell2mat(databank.data(:,(strcmp(databank.nm,'FeatureSample'))==1));
    StimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd = cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    StimStrength = StimStrength(isnan(data)==0);
    pstd = pstd(isnan(data)==0);
end

%----------------------------------
%1 - Sort data as training and test
%----------------------------------
%Divide data into 5 sets, gets best fit parameters on 4 sets (training) and
%use them to calculate the R2 of the model fitted to the test data.
fprintf('%s \n','(CrossValR2fit) Sorting training and data sets')
numSet = 5;
output.CrossVR2 = nan(numSet,1);
%loop over sets
for thisSet = 1 : numSet
    
    %status
    fprintf('%s \n',['(CrossValR2fit) Calculating R-squared for set ',...
        num2str(thisSet),'/',num2str(numSet)])
    
    %------------------------------------------------------------
    %1 - Divide data into 5 sets that each contain all conditions
    %------------------------------------------------------------
    %(case von Mises)
    %----------------
    if strcmp(priorShape,'vonMisesPrior')
        
        %each trial's experimental condition
        [~,~,idxCondAlltrials] = SLuniqpair([pstd StimStrength disp]);
        numC = max(idxCondAlltrials);
        output.CrossVSet = nan(size(data));
        output.CrossVnumSet = numSet;
        setID = 1 : output.CrossVnumSet;
        
        for Ci = 1 : max(numC)
            
            %separate each single condition in 5 sets of n trials.
            TrialsThisC=Ci==idxCondAlltrials;
            numTrialsThisC=sum(TrialsThisC);
            numTrialsPerSetThisC=fix(numTrialsThisC/numSet);
            
            %assign a set to each trial
            sets=setID(ones(numTrialsPerSetThisC,1),:);
            sets=sets(:);
            sets(end:numTrialsThisC)=numSet;
            output.CrossVSet(TrialsThisC)=sets;
        end
        
        %training trials (4 sets, 5th is left out for test)
        testSetID=thisSet;
        TrialsTrain=output.CrossVSet~=testSetID;
        
        %training data and conditions
        dataTrain=data(TrialsTrain);
        dispTrain=disp(TrialsTrain);
        StimStrengthTrain=StimStrength(TrialsTrain);
        pstdTrain=pstd(TrialsTrain);
        priorModesTrain=priorModes(TrialsTrain);
    end
    
    %------------------------------
    %2 - fit model to training data
    %------------------------------
    %(case von Mises)
    %----------------
    %Calculate data mean and std
    if strcmp(priorShape,'vonMisesPrior')
        
        %status
        fprintf('%s \n',['(CrossValR2fit) Calculating training data means'...
            ' and stds per condition for set ',...
            num2str(thisSet),'/',num2str(numSet)])
        
        %data are sorted by experimental conditions
        [meanData,stdData,myF] = SLCircStat(dataTrain,dispTrain,StimStrengthTrain,pstdTrain);
        meanData=meanData(:);
        stdData=stdData(:);
        myF1=myF.f1.D3(:);
        myF2=myF.f2.D3(:);
        myF3=myF.f3.D3(:);
        myCond=[myF3 myF2 myF1];
        
        %data descriptive stats
        statsTrain = SLmakeCircStat(dataTrain,pstdTrain,StimStrengthTrain,dispTrain);
        
        %status
        fprintf('%s \n',['(CrossValR2fit) IMPORTANT !!!! Please Check',...
            ' the training data that is plotted. e.g.,number of sample data for',...
            ' a condition may not be enough. Not enough data will provide',...
            ' poor R-squared ! The reason is that there may be a single',...
            ' trial at priors tails that produce 0 deg std. This is wrong !',...
            ' make sure there are at least 3 trials at the tails'])
        
        %display descriptive stats to check before fitting
        %M = [statsTrain.count statsTrain.conditions];
        %fprintf('%8.0f %8.0f %8.2f %8.0f \n',M')
        
        %back up
        output.statsTrain{thisSet} = statsTrain;
        
        %ckeck data mean and std visually
        %fprintf('%s \n',['(CrossValR2fit) IMPORTANT !!!! Now drawing data',...
        %    ' from training set....'])
        %drawMeanPre(meanData,stdData,[],[],myCond);
    end
    
    %for fminsearch - Nelder-Mead
    %OptionsNM=optimset('Display','iter');
    OptionsNM = [];
    
    %fit if no input parameters
    if isempty(fitPinput)==1
        
        %Sets of initial parameters: best matching parameters are found by
        %matching raw data with simulations graphically (called "Kreal").
        %with best matching initial values
        kllh = initp(1:3);
        klearnt = initp(4:7);
        kcardinal = initp(8);
        fractRand = initp(9);
        motorNoise = initp(10);
        
        %(1) best matching(true) parameters
        k0(1,:) = [kllh  klearnt  kcardinal fractRand motorNoise];
        %(2) 8x stronger priors & best matching llh
        k0(2,:) = [kllh    8*klearnt    kcardinal fractRand motorNoise ];
        %(3) 8x weaker priors & best matching llh
        k0(3,:) = [kllh    klearnt./8   kcardinal fractRand motorNoise ];
        %(4) 8x stronger likelihoods & best matching priors
        k0(4,:) = [kllh.*8 klearnt      kcardinal fractRand motorNoise ];
        %(5) 8x stronger likelihoods & stronger priors
        k0(5,:) = [[kllh   klearnt].*8  kcardinal fractRand motorNoise ];
        %(6) 8x stronger likelihoods & weaker priors
        k0(6,:) = [kllh.*8 klearnt./8   kcardinal fractRand motorNoise ];
        %(7) 8x weaker likelihood & true priors
        k0(7,:) = [kllh./8 klearnt      kcardinal fractRand motorNoise ];
        %(8) 8x weaker likelihood & stronger priors
        k0(8,:) = [kllh./8 klearnt.*8   kcardinal fractRand motorNoise ];
        %(9) 8x weaker likelihood & weaker priors
        k0(9,:) = [[kllh   klearnt]./8  kcardinal fractRand motorNoise ];
        %(10)likelihoods & priors are same.
        k0(10,:) =[repmat(nanmean(k0(1,1:8)),1,7) kcardinal fractRand motorNoise];
        
        %status
        output.numInitP = size(k0,1);
        fprintf('%s \n',['(CrossValR2fit) Initialize your',...
            num2str(output.numInitP),' initial parameters'])
        
        %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            TheModel = 'withoutCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2fit) The model does not contain',...
                'a cardinal prior'])
        else
            TheModel = 'withCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2fit) The model contains',...
                'a cardinal prior'])
        end
        
        %fitting
        %k can be anything >=0;
        fitPbkp=nan(size(k0,1),size(k0,2));
        SSE_bkp=nan(size(k0,1),1);
        Hess = [];
        
        %options
        output.numfitP = size(k0,2);
        output.initParameters = k0;
        
        %for default fitting
        %options = [];
        
        %for fast debugging or fitting
        output.numIter = 100;
        output.numFunEval = 100;
        options = optimset('MaxIter',output.numIter,'MaxFunEvals',output.numFunEval);
        fprintf('%s \n',['(CrossValR2fit) Fitting with ',num2str(output.numIter),...
            ' iterations and ',num2str(output.numFunEval),' function evaluations...'])
        vrg = varargin{1};
        
        parfor i=1:size(k0,1)
            %for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n',['(CrossValR2fit)',' Set of initial parmeters: ',...
                num2str(i),'/',num2str(size(k0,1))])
            
            %Nelder-Mead
            [fitPtmp,SSE,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                makeSSE(dispTrain,...
                StimStrengthTrain,...
                pstdTrain,...
                fitPtmp,...
                priorShape,...
                priorModesTrain,...
                TheModel,...
                meanData,stdData,myCond,...
                vrg),...
                k0(i,:),...
                options);
            
            %Fit parameters and SSE
            fitPbkp(i,:) = fitPtmp;
            SSE_bkp(i) = SSE;
            outputFitBkp(i) = {outputFit};
            exitflagBkp{i}  = exitflag;
        end
        
        %fit backups
        output.fitPbkp{thisSet}     = fitPbkp;
        output.SSE_bkp{thisSet}     = SSE_bkp;
        output.outputFit{thisSet}   = outputFitBkp;
        output.fitexitflag{thisSet} = exitflagBkp;
        
        %Max likelihood. To have an idea of what log likelihood we should
        %obtain our reasoning is that if one cannot predict at all motion
        %direction at a given trial, the probability that any given direction
        %be chosen is 1/360. (we consider estimate with 1 degree resolution).
        %So loglikelihood should be log(1/360)=-5.88. This mutiplied (because
        %of log) by our ~5500 trials, we get negLogl~32000. Any ability to
        %somewhat predict better motion direction better than chance yields
        %soemthing lower than 32000.
        fprintf('%s \n',['(CrossValR2fit) Identifying minimum SSE and',...
            ' best fit parameters'])
        [minSSE,position] = min(SSE_bkp(:));
        i = ind2sub(size(SSE_bkp),position);
        
        %minimal SSE
        output.minSSE = minSSE;
        
        %optimal parameters
        output.fitP.p(thisSet,:) = fitPbkp(i,:);
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'vonMisesPrior')
            output.fitP.nm={'StimStrength24','StimStrength12','StimStrength6','pstd80','pstd40','pstd20','pstd10',...
                'kcardinal','Prand','km'};
        end
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'bimodalPrior')
            output.fitP.nm={'StimStrength24','StimStrength12','StimStrength6','pstd145_305','pstd165_285',...
                'pstd185_265','pstd205_245',...
                'kcardinal','Prand','km'};
        end
        
        %(case we did not fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            output.fitP.p(thisSet,strcmp(output.fitP.nm,'kcardinal'))=NaN;
        end
        
        %predictions
        pred=[];
        
        %std of model parameters for interior-point only.
        stdPa=[];
        
        %save data in case bug
        save(filename)
        
        %if model fit parameters are input manually just give log likelihood of
        %data, and make predictions based on those fit parameters.
    elseif isempty(fitPinput)==0
        
        %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
        if isnan(fitPinput(8))
            TheModel='withoutCardinal';
        else
            TheModel='withCardinal';
        end
        
        SSE = makeSSE(dispTrain,StimStrengthTrain,pstdTrain,fitPinput,...
            priorShape,priorModesTrain,TheModel,meanData,stdData,...
            dataCond,varargin{1});
        
        %SSE
        output.minSSE = SSE;
        
        %fitP
        output.fitPbkp{thisSet} = fitPinput;
        output.fitP.p(thisSet,:) = fitPinput;
        output.fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','klearnt80','klearnt40','klearnt20',...
            'klearnt10','Kcardinal','Prand','km'};
        pred = [];
        stdPa = [];
    end
    
    %----------------------------------------------------------------------
    %2 -  Use training best fit parameters to calculate R2 in test data set
    %----------------------------------------------------------------------
    
    %status
    fprintf('%s \n',['(CrossValR2fit) Use training best fit parameters to',...
        'Calculate R-squared in test data...'])
    
    %make test data set
    %select test trials (5th set left out for test)
    output.TrialsTest{thisSet} = output.CrossVSet==testSetID;
    dataTest = data(output.TrialsTest{thisSet});
    
    %test conditions
    dispTest = disp(output.TrialsTest{thisSet});
    StimStrengthTest  = StimStrength(output.TrialsTest{thisSet});
    pstdTest = pstd(output.TrialsTest{thisSet});
    priorModesTest = priorModes(output.TrialsTest{thisSet});
    
    %mean and std of test data
    [meanDataTest,stdDataTest,myFTest]=SLCircStat(dataTest,dispTest,...
        StimStrengthTest,pstdTest);
    meanDataTest = meanDataTest(:);
    stdDataTest = stdDataTest(:);
    myF1 = myFTest.f1.D3(:);
    myF2 = myFTest.f2.D3(:);
    myF3 = myFTest.f3.D3(:);
    myCondTest = [myF3 myF2 myF1];
    
    
    %data descriptive stats
    statsTest = SLmakeCircStat(dataTest,pstdTest,StimStrengthTest,dispTest);
    
    %status
    fprintf('%s \n',['(CrossValR2fit) IMPORTANT !!!! Please Check',...
        ' the test data that is plotted. e.g.,number of sample data for',...
        ' a condition may not be enough. Not enough data will provide',...
        ' poor R-squared !'])
    
    %display descriptive stats to check before fitting
    %M = [statsTest.count statsTest.conditions];
    %fprintf('%8.0f %8.0f %8.2f %8.0f \n',M')
    
    %back up
    output.statsTest{thisSet} = statsTest;
    
    %ckeck data mean and std visually
    fprintf('%s \n',['(CrossValR2fit) IMPORTANT !!!! Now drawing data',...
        ' from test set....'])
    SLdrawModelsPredictionCentered(meanDataTest,stdDataTest,[],[],myCondTest,...
        priorModes,priorShape,'yCentered')
    
    %test SSE
    SSE = makeSSE(dispTest,StimStrengthTest,pstdTest,output.fitP.p(thisSet,:),priorShape,...
        priorModesTest,TheModel,meanDataTest,stdDataTest,myCondTest,varargin{1});
    
    %test R-squared
    output.CrossVR2(thisSet) = makeR2(meanDataTest,stdDataTest,SSE);
    
    %status
    fprintf('%s \n','----------------------------------------------------')
    fprintf('%s \n',['(CrossValR2fit) R-squared for set ',...
        num2str(thisSet),'/',num2str(numSet),'= ',...
        num2str(output.CrossVR2(thisSet))])
    fprintf('%s \n','----------------------------------------------------')
end

%model's mean R2
%status
fprintf('%s \n',['(CrossValR2fit) Calculating average R2 over the',...
    num2str(numSet),' test sets'])

output.CrossVmeanR2 = nanmean(output.CrossVR2);

%cross-validated R^2 fit (trial-data)
function [pred,data,disp,StimStrength,pstd,stdPa,output] = CrossValR2FitTrialData(...
    databank,...
    fitPinput,...
    initp,...
    priorShape,...
    filename,...
    output,...
    varargin)

%STEPS of the procedure
%1. First, sort data as training and test.

%2. Second, get the best fit parameters that minimize sum of square error
%between training data and predictions.

%3. Use the fit Parameters to generate predictions for the test data set
%and calculate test R2.

%4. Start a new training data set and repeat the procedure.

%5. Store the R2 for each of the 5 sets and calculate the mean R2.

%status
fprintf('%s \n','(CrossValR2FitTrialData) Now fitting model to trial data')

%(case check fitting fake data).
%Get simulated data in the workspace.
if sum(strcmp(varargin{1},'testFakedata'))
    data = evalin('base','pred');
    disp = evalin('base','d');
    StimStrength = evalin('base','StimStrength');
    pstd = evalin('base','pstd');
    priorModes = evalin('base','priorModes');
    
else
    
    %subjects' data
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'estimatedFeature'))==1)));
    
    %conditions
    disp = cell2mat(databank.data(:,(strcmp(databank.nm,'FeatureSample'))==1));
    StimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd  =cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    StimStrength = StimStrength(isnan(data)==0);
    pstd = pstd(isnan(data)==0);
end

%store data
output.data = data;
output.disp = disp;
output.pstd = pstd;
output.priorModes = priorModes;

%----------------------------------
%1 - Sort data as training and test
%----------------------------------
%Divide data into 2 sets, gets best fit parameters on 1 sets (training) and
%use them to calculate the R2 of the model fitted to the test data.
fprintf('%s \n','(CrossValR2FitTrialData) Sorting training and test sets ...')
numSet = 2;
output.CrossVR2 = nan(numSet,1);

%loop over sets
for thisSet = 1 : numSet
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialData) Calculating R^2 for set ',...
        num2str(thisSet),'/',num2str(numSet)])
    
    %------------------------------------------------------------
    %1 - Divide data into 2 sets that each contain all conditions
    %------------------------------------------------------------
    %(case von Mises)
    %----------------
    if strcmp(priorShape,'vonMisesPrior')
        
        %each trial's experimental condition
        [~,~,idxCondAlltrials] = SLuniqpair([pstd StimStrength disp]);
        numC = max(idxCondAlltrials);
        output.CrossVSet = nan(size(data));
        output.CrossVnumSet = numSet;
        setID = 1:output.CrossVnumSet;
        
        for Ci = 1 : max(numC)
            
            %separate each single condition in 5 sets of n trials.
            TrialsThisC = Ci==idxCondAlltrials;
            numTrialsThisC = sum(TrialsThisC);
            numTrialsPerSetThisC = fix(numTrialsThisC/numSet);
            
            %assign a set to each trial
            sets = setID(ones(numTrialsPerSetThisC,1),:);
            sets = sets(:);
            sets(end:numTrialsThisC) = numSet;
            output.CrossVSet(TrialsThisC) = sets;
        end
        
        %training trials (1 set, 2nd is left out for test)
        testSetID = thisSet;
        TrialsTrain = output.CrossVSet~=testSetID;
        
        %training data and conditions
        dataTrain = data(TrialsTrain);
        dispTrain = disp(TrialsTrain);
        StimStrengthTrain = StimStrength(TrialsTrain);
        pstdTrain = pstd(TrialsTrain);
        priorModesTrain = priorModes(TrialsTrain);
    end
    
    %store
    output.dataTrainBkp{thisSet} = dataTrain;
    
    %------------------------------
    %2 - fit model to training data
    %------------------------------
    %(case von Mises)
    %----------------
    %Calculate data mean and std
    if strcmp(priorShape,'vonMisesPrior')
        
        %status
        fprintf('%s \n',['(CrossValR2FitTrialData) Calculating training data means'...
            ' and stds per condition for set ',...
            num2str(thisSet),'/',num2str(numSet)])
        
        %training data mean and std are sorted by conditions
        [meanData,stdData,myF] = SLCircStat(dataTrain,dispTrain,...
            StimStrengthTrain,pstdTrain);
        meanData = meanData(:);
        stdData = stdData(:);
        myF1 = myF.f1.D3(:);
        myF2 = myF.f2.D3(:);
        myF3 = myF.f3.D3(:);
        myCond = [myF3 myF2 myF1];
        
        %data descriptive stats
        %statsTrain = SLmakeCircStat(dataTrain,pstdTrain,StimStrengthTrain,dispTrain);
        
        %status
        %fprintf('%s \n',['(CrossValR2FitTrialData) IMPORTANT !!!! Please Check',...
        %   ' the training data that is plotted. e.g.,number of sample data for',...
        %   ' a condition may not be enough. Not enough data will provide',...
        %   ' poor R-squared ! The reason is that there may be a single',...
        %   ' trial at priors tails that produce 0 deg std. This is wrong !',...
        %   ' make sure there are at least 3 trials at the tails'])
        
        %display descriptive stats to check before fitting
        %M = [statsTrain.count statsTrain.conditions];
        %fprintf('%8.0f %8.0f %8.2f %8.0f \n',M')
        
        %back up
        %output.statsTrain{thisSet} = statsTrain;
        
        %ckeck data mean and std visually
        %fprintf('%s \n',['(CrossValR2FitTrialData) IMPORTANT !!!! Now drawing data',...
        %    ' from training set....'])
        %drawMeanPre(meanData,stdData,[],[],myCond);
    end
    
    %for fminsearch - Nelder-Mead
    %OptionsNM=optimset('Display','iter');
    OptionsNM = [];
    
    %fit if no input parameters
    if isempty(fitPinput)==1
        
        %Sets of initial parameters: best matching parameters are found by
        %matching raw data with simulations graphically (called "Kreal").
        %We pre-generate 27 reasonable sets of initial parameters.
        %+-----+-------+----------------------+
        %|     |       |          prior       |
        %.-----|-------+-------+-------+------+
        %|     |       | true  | strong| weak |
        %.-----|-------+-------+-------+------+
        %|     |true   |(1)t-t |(2)t-s |(3)t-w|
        %.     |-------+-------+-------+------+
        %| llh |strong |(4)s-t |(5)s-s |(6)s-w|
        %|     |-------+-------+-------+------+
        %|     |weak   |(7)w-t |(8)w-s |(9)w-w|
        %+-----+-----------------------+------+
        
        %set 10: likelihood and prior strengths are all same.
        %eventually we may use later
        %10 sets.
        
        %note: strong priors (and weak priors) are 8 times stronger than best
        %matching priors. 8x is the factor for which I see clear deviation
        %of simulation from data.
        
        %We used one intial value for probability of random estimation, motor
        %noise and cardinal prior strength, that we think are relatively small
        %values (high k for motor noise is low motor noise) by looking at the
        %data.
        
        %best graphically matching initial values found by manual
        %fitting
        kllh = initp(1:3);
        klearnt = initp(4:7);
        kcardinal = initp(8);
        fractRand = initp(9);
        motorNoise = initp(10);
        
        %(1) best matching(true) parameters
        k0(1,:)=[kllh    klearnt      kcardinal fractRand motorNoise];
        %(2) 8x stronger priors & best matching llh
        k0(2,:)=[kllh    8*klearnt    kcardinal fractRand motorNoise ];
        %(3) 8x weaker priors & best matching llh
        k0(3,:)=[kllh    klearnt./8   kcardinal fractRand motorNoise ];
        %(4) 8x stronger likelihoods & best matching priors
        k0(4,:)=[kllh.*8 klearnt      kcardinal fractRand motorNoise ];
        %(5) 8x stronger likelihoods & stronger priors
        k0(5,:)=[[kllh   klearnt].*8  kcardinal fractRand motorNoise ];
        %(6) 8x stronger likelihoods & weaker priors
        k0(6,:)=[kllh.*8 klearnt./8   kcardinal fractRand motorNoise ];
        %(7) 8x weaker likelihood & true priors
        k0(7,:)=[kllh./8 klearnt      kcardinal fractRand motorNoise ];
        %(8) 8x weaker likelihood & stronger priors
        k0(8,:)=[kllh./8 klearnt.*8   kcardinal fractRand motorNoise ];
        %(9) 8x weaker likelihood & weaker priors
        k0(9,:)=[[kllh   klearnt]./8  kcardinal fractRand motorNoise ];
        %(10)likelihoods & priors are same.
        k0(10,:)=[repmat(nanmean(k0(1,1:8)),1,7) kcardinal fractRand motorNoise];
        
        %status
        output.numInitP = size(k0,1);
        fprintf('%s \n',['(CrossValR2FitTrialData) Initialize your ',...
            num2str(output.numInitP),' initial parameters'])
        
        %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            TheModel = 'withoutCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialData) The model does not contain',...
                'a cardinal prior'])
        else
            TheModel = 'withCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialData) The model contains',...
                'a cardinal prior'])
        end
        
        %fitting
        %k can be anything >=0;
        fitPbkp = nan(size(k0,1),size(k0,2));
        SSE_bkp = nan(size(k0,1),1);
        Hess = [];
        
        %options
        output.numfitP = size(k0,2);
        output.initParameters = k0;
        
        %for default fitting
        %options = [];
        
        %for fast debugging or fitting
        output.numIter = 100;
        output.numFunEval = 100;
        options = optimset('MaxIter',output.numIter,'MaxFunEvals',...
            output.numFunEval);
        
        %status
        fprintf('%s \n',['(CrossValR2FitTrialData) Fitting with ',...
            num2str(output.numIter),...
            ' iterations and ',num2str(output.numFunEval),...
            ' function evaluations...'])
        
        %loop over 10 sets of initial parameters
        vrg = varargin{1};
        
        parfor i = 1 : size(k0,1)
            %for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialData)',' Set of initial parmeters: ',...
                num2str(i),'/',num2str(size(k0,1))])
            
            %Nelder-Mead (simplex search algorithm)
            [fitPtmp,SSE,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                makeSSE(dispTrain,...
                StimStrengthTrain,...
                pstdTrain,...
                fitPtmp,...
                priorShape,...
                priorModesTrain,...
                TheModel,...
                'Trial',...
                dataTrain,...
                meanData,stdData,myCond,...
                output,...
                vrg),...
                k0(i,:),...
                options);
            
            %Fit parameters and SSE
            fitPbkp(i,:)    = fitPtmp;
            SSE_bkp(i)      = SSE;
            outputFitBkp(i) = {outputFit};
            exitflagBkp{i}  = exitflag;
        end
        
        %backups
        output.fitPbkp{thisSet}     = fitPbkp;
        output.SSE_bkp{thisSet}     = SSE_bkp;
        output.outputFit{thisSet}   = outputFitBkp;
        output.fitexitflag{thisSet} = exitflagBkp;
        
        %Max likelihood. To have an idea of what log likelihood we should
        %obtain our reasoning is that if one cannot predict at all motion
        %direction at a given trial, the probability that any given direction
        %be chosen is 1/360. (we consider estimate with 1 degree resolution).
        %So loglikelihood should be log(1/360)=-5.88. This mutiplied (because
        %of log) by our ~5500 trials, we get logL~32000. Any ability to
        %somewhat predict better motion direction better than chance yields
        %soemthing lower than 32000.
        %status
        %status
        fprintf('%s \n',['(CrossValR2FitTrialData) Identifying minimum SSE and',...
            ' best fit parameters'])
        
        [minSSE,position] = min(SSE_bkp(:));
        i = ind2sub(size(SSE_bkp),position);
        
        %minimal SSE
        output.minSSE=minSSE;
        
        %optimal parameters
        output.fitP.p(thisSet,:)=fitPbkp(i,:);
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'vonMisesPrior')
            output.fitP.nm={'StimStrength24','StimStrength12','StimStrength6','pstd80','pstd40','pstd20','pstd10',...
                'kcardinal','Prand','km'};
        end
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'bimodalPrior')
            output.fitP.nm={'StimStrength24','StimStrength12','StimStrength6','pstd145_305','pstd165_285',...
                'pstd185_265','pstd205_245',...
                'kcardinal','Prand','km'};
        end
        
        %(case we did not fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            output.fitP.p(thisSet,strcmp(output.fitP.nm,'kcardinal'))=NaN;
        end
        
        %predictions
        pred=[];
        
        %std of model parameters for interior-point only.
        stdPa=[];
        
        %save data in case bug
        save(filename)
        
        %if model fit parameters are input manually just give log likelihood of
        %data, and make predictions based on those fit parameters.
    elseif isempty(fitPinput)==0
        
        %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
        if isnan(fitPinput(8))
            TheModel='withoutCardinal';
        else
            TheModel='withCardinal';
        end
        
        %SSE
        SSE = makeSSE(dispTrain,StimStrengthTrain,pstdTrain,fitPinput,...
            priorShape,priorModesTrain,TheModel,'Trial',dataTrain,meanData,...
            stdData,dataCond,output,varargin{1});
        
        %SSE
        output.minSSE = SSE;
        
        %fitP
        output.fitPbkp{thisSet} = fitPinput;
        output.fitP.p(thisSet,:) = fitPinput;
        output.fitP.nm = {'StimStrength24','StimStrength12','StimStrength6','klearnt80','klearnt40','klearnt20',...
            'klearnt10','Kcardinal','Prand','km'};
        pred=[];
        stdPa=[];
    end
    
    %----------------------------------------------------------------------
    %2 -  Use training best fit parameters to calculate R2 in test data set
    %----------------------------------------------------------------------
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialData) Use training best fit parameters to',...
        'Calculate R-squared in test data...'])
    
    %make test data set
    %select test trials (5th set left out for test)
    output.TrialsTest{thisSet} = output.CrossVSet==testSetID;
    dataTest = data(output.TrialsTest{thisSet});
    output.dataTestBkp = dataTest;
    
    %test conditions
    dispTest = disp(output.TrialsTest{thisSet});
    StimStrengthTest  = StimStrength(output.TrialsTest{thisSet});
    pstdTest = pstd(output.TrialsTest{thisSet});
    priorModesTest = priorModes(output.TrialsTest{thisSet});
    
    %mean and std of test data
    [meanDataTest,stdDataTest,myFTest] = SLCircStat(dataTest,dispTest,...
        StimStrengthTest,pstdTest);
    meanDataTest = meanDataTest(:);
    stdDataTest = stdDataTest(:);
    myF1 = myFTest.f1.D3(:);
    myF2 = myFTest.f2.D3(:);
    myF3 = myFTest.f3.D3(:);
    myCondTest = [myF3 myF2 myF1];
    SLdrawModelsPredictionCentered(meanDataTest,stdDataTest,[],[],myCondTest,...
        priorModes,priorShape,'yCentered')
    
    
    %data descriptive stats
    statsTest = SLmakeCircStat(dataTest,pstdTest,StimStrengthTest,dispTest);
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialData) IMPORTANT !!!! Please Check',...
        ' the test data that is plotted. e.g.,number of sample data for',...
        ' a condition may not be enough. Not enough data will provide',...
        ' poor R-squared !'])
    
    %display descriptive stats to check before fitting
    %M = [statsTest.count statsTest.conditions];
    %fprintf('%8.0f %8.0f %8.2f %8.0f \n',M')
    
    %back up
    output.statsTest{thisSet} = statsTest;
    
    %ckeck data mean and std visually
    fprintf('%s \n',['(CrossValR2FitTrialData) IMPORTANT !!!! Now drawing data',...
        ' from test set....'])
    SLdrawModelsPredictionCentered(meanDataTest,stdDataTest,[],[],myCondTest,...
        priorModes,priorShape,'yCentered')
    
    %test SSE
    SSE = makeSSE(dispTest,StimStrengthTest,pstdTest,output.fitP.p(thisSet,:),priorShape,...
        priorModesTest,TheModel,'Trial',dataTest,meanDataTest,stdDataTest,...
        myCondTest,output,varargin{1});
    
    %test R-squared
    output.CrossVR2(thisSet) = makeR2('Trial',meanDataTest,stdDataTest,...
        dataTest,SSE);
    
    %status
    fprintf('%s \n','----------------------------------------------------')
    fprintf('%s \n',['(CrossValR2FitTrialData) R-squared for set ',...
        num2str(thisSet),'/',num2str(numSet),'= ',...
        num2str(output.CrossVR2(thisSet))])
    fprintf('%s \n','----------------------------------------------------')
end

%mean R^2
%status
fprintf('%s \n',['(CrossValR2FitTrialData) Calculating average R2 over the',...
    num2str(numSet),' test sets'])

output.CrossVmeanR2 = mean(output.CrossVR2);


%SSE between model & data (trial & average)
function [SSE,fitP,output] = makeSSE(displ,StimStrength,pstd,fitP,priorShape,priorModes,...
    TheModel,TrialOrMean,trialData,meanData,stdData,dataCond,output,varargin)


%time
ticSSE = tic;

%remove missing conditions (NaN) in the data
pos = ~isnan(meanData);
meanData = meanData(pos);
stdData = stdData(pos);
dataCond = dataCond(pos,:);

%(case cardinal priors)
if isnan(fitP(8))
    fitP(8) = 0;
end

%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time. Constrain other fraction
%parameters between 0 and 1, and other parameters as>0. Case a parameter
%that is not cardinal prior strength is missing (NaN) fraction of random
%estimates
if fitP(9) > 1
    SSE = 1e9;
    return
end
if any(fitP < 0)
    SSE = 1e9;
    return
end
if any(isnan(fitP([1,2,3,4,5,6,7,9,10])))
    SSE = 1e9;
    fprintf('%s \n',['(makeSSE) One of your fit parameter',...
        'that is not Kcardinal is NaN'])
    keyboard
end

%predictions
[meanPred,stdPred,cond,~,~,output] = SLmakePredictionsBayesianModel([],displ,StimStrength,pstd,fitP,...
    priorShape,priorModes,TheModel,TrialOrMean,output,varargin{1});



%case least-square fit to estimate mean and std
%----------------------------------------------
if strcmp(TrialOrMean,'Mean')
    
    %check if data and predictions are ordered in the same way
    if sum(sum(dataCond - Predcond)) ~= 0
        fprintf('%s \n',['(makeSSE) Something s wrong. Data and predictions are not ',...
            'ordered in the same way'])
        keyboard
    end
    
    %SSE between model predictions and data
    %circular error
    ErrorForMean = SLvectors2signedAngle(meanData,meanPred);
    
    %linear error
    ErrorForstd = stdData - stdPred;
    
    %error
    Error = [ErrorForMean; ErrorForstd];
    
    %%(Visual check) estimate mean and std and predictions
    %SLdrawModelsPredictionCentered(meanData,stdData,meanPred,stdPred,dataCond,'yCentered')
end

%case least-square fit to trial data
%-----------------------------------
if strcmp(TrialOrMean,'Trial')
    trialPred = output.TrialPred;
    Error = SLvectors2signedAngle(trialData,trialPred);
    
    %%(Visual check) data mean and pred
    %%data
    %[meanData,stdData,myF] = SLCircStat(trialData,displ,...
    %   StimStrength,pstd);
    %meanData = meanData(:);
    %stdData = stdData(:);
    %myCond = [myF.f3.D3(:) myF.f2.D3(:) myF.f1.D3(:)];
    
    %%predictions
    %%predictions of std are not possible because model outputs the same mean
    %%for each same condition.
    %meanPred = SLCircStat(trialPred,displ,...
    %   StimStrength,pstd);
    %meanPred = meanPred(:);
    
    %%draw
    %SLdrawModelsPredictionCentered(meanData,stdData,meanPred,[],myCond,'yCentered')
end

%SSE
SSE = nansum(Error.^2);

%stop to debug if SSE is a complex value
if ~isreal(SSE)
    dbstop
end

%print fitting info
ti = toc(ticSSE);
fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
    SSE,fitP,ti)

%R^2 (trial & average)
function R2 = makeR2(TrialOrMean,meanData,stdData,trialData,SSE)
%R-squared = 1 - SSE/SST
%SSE (alias residual) is the SSE calculated for the best model's prediction.
%note: R^ (alias coefficient of determination) can be <0 when the model
%does an awful job predicting the data. In this case SSE exceeds that is
%the model fits the data even worse than does a horizontal line.
%Model may not be appropriate or constraints may not be set correctly.
%see http://www.graphpad.com/support/faqid/711/
%
%Another way is:
%function R2=makeR22(data,pred)
%%R2=1 - SSE/SST
%%SSE (alias residual) is the SSE calculated for the best model's prediction.
%%note: R^ (alias coefficient of determination) can be <0 when the model
%%does an awful job predicting the data. In this case SSE exceeds that is
%%the model fits the data even worse than does a horizontal line.
%%Model may not be appropriate or constraints may not be set correctly.
%%see http://www.graphpad.com/support/faqid/711/
%
%R=corr(data,pred);
%R2=R^2;

%case estimate mean and std
%--------------------------
if strcmp(TrialOrMean,'Mean')
    
    %circular mean data mean and mean data std
    stat1 = SLcircMeanStd(meanData);
    MeanOfMeanData = stat1.deg.mean;
    stat2 = SLcircMeanStd(stdData);
    MeanOfStdData = stat2.deg.mean;
    
    %errors
    errorForMean = SLvectors2signedAngle(meanData,MeanOfMeanData);
    errorForStd = SLvectors2signedAngle(stdData,MeanOfStdData);
    error = [errorForMean; errorForStd];
end

%case trial data
%---------------
if strcmp(TrialOrMean,'Trial')
    
    %circular data mean
    stat = SLcircMeanStd(trialData);
    meanDataAll = stat.deg.mean;
    
    %circular error
    error = SLvectors2signedAngle(trialData,meanDataAll);
end

%SST
SST = nansum(error.^2);

%R^2
R2 = 1 - SSE/SST;


%-------------------
%Nested calculations
%-------------------
%data mean,std per condition
function [meanData,stdData,dataCond] = makeDataMeanAndStd(data,disp,...
    StimStrength,pstd,priormodes,priorShape)

%(case von Mises prior)
%----------------------
%Calculate data mean and std for each experimental condition (disp,StimStrength,pstd)
if strcmp(priorShape,'vonMisesPrior')
    
    %data are sorted by experimental conditions
    [meanData,stdData,myF] = SLCircStat(data,disp,StimStrength,pstd);
    meanData = meanData(:);
    stdData = stdData(:);
    myF1 = myF.f1.D3(:);
    myF2 = myF.f2.D3(:);
    myF3 = myF.f3.D3(:);
    dataCond = [myF3 myF2 myF1];
end

%(case bimodal prior)
%----------------------
if strcmp(priorShape,'bimodalPrior')
    
    %prior conditions
    priorCond = priormodes(:,2) - priormodes(:,1);
    PriorCondunq = unique(priorCond);
    
    %check that we get all prior conditions
    numPriorCond = size(SLuniqpair(priormodes),1);
    if numel(PriorCondunq) == numPriorCond
        
        %data are sorted by experimental conditions
        [meanData,stdData,myF] = SLCircStat(data,disp,StimStrength,priorCond);
        meanData = meanData(:);
        stdData = stdData(:);
        myF1 = myF.f1.D3(:);
        myF2 = myF.f2.D3(:);
        myF3 = myF.f3.D3(:);
        dataCond = [myF3 myF2 myF1];
        
        %if not
    else
        fprintf('%s \n',['(makeDataMeanAndStd) Something wrong with number',...
            ' of prior conditions'])
    end
    
end

%remove missing conditions (NaN) in the data
pos =~ isnan(meanData);
meanData = meanData(pos);
stdData = stdData(pos);
dataCond = dataCond(pos,:);

%data distributions per condition
function [p,xpdf] = makeDataDist(data,d,StimStrength,pstd,priorModes,cond,priorShape)

%need to start from 0 and end at 360 to get the full range
motDir = 0:10:360;

%this stimStrength
p = nan(numel(motDir)-1,size(cond,1));

%prior conditions
%(case bimodal prior)
%--------------------
if strcmp(priorShape,'bimodalPrior')
    
    priorCond = priorModes(:,2) - priorModes(:,1);
    
    %sanity check
    numPriorCond = size(SLuniqpair(priorModes),1);
    if numel(unique(priorCond)) ~= numPriorCond
        fprintf('%s \n', ['(makeDataDist) Something wrong with number of',...
            'prior conditions for bimodal prior'])
    end
end

%this condition
YesNo = input('Do you want to see sample size for each condition y/n ?','s');

for i = 1 : size(cond,1)
    
    %this condition
    %case von Mises prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        thisCon = pstd==cond(i,1) & StimStrength==cond(i,2) & d==cond(i,3);
        
        %case bimodal prior
        %--------------------
    elseif strcmp(priorShape,'bimodalPrior')
        thisCon = priorCond==cond(i,1) & StimStrength==cond(i,2) & d==cond(i,3);
        
    else
        fprintf('%s \n', '(makeDataDist) You need to input prior type')
    end
    
    %get data for this condition
    dtoHist = data(thisCon);
    
    %make probability distribution
    [~,xpdf,ctmp(:,i)] = makePdf(motDir,dtoHist,'raw');
    
    %adjust count
    c(:,i) = ctmp(1:end-1,i);
    
    %probability
    p(:,i) = c(:,i)/sum(c(:,i));
    %[p(:,i),xpdf,c(:,i)] = makePdf(motDir,dtoHist,'raw');
    
    %if show sample size
    if strcmp(YesNo,'y')
        fprintf('%s \n', ['(makeDataDist) ',num2str(numel(dtoHist)),' trials'])
    end
    
end


%-----------------
%drawing functions
%-----------------
%model mean, std
function drawMeanPre(meanData,stdData,meanPred,stdPred,dataCond)

%make sure column vector
if size(meanPred,2) > size(meanPred,1)
    meanPred=meanPred';
end
if size(meanData,2) > size(meanData,1)
    meanData=meanData';
end
if size(stdData,2) > size(stdData,1)
    stdData=stdData';
end
if size(stdPred,2) > size(stdPred,1)
    stdPred=stdPred';
end

%case meanData and stdData are not input
if isempty(meanData)==1
    meanData = nan(size(dataCond,1),1);
end
if isempty(stdData)==1
    stdData = nan(size(dataCond,1),1);
end
if isempty(meanPred)==1
    meanPred = nan(size(dataCond,1),1);
end
if isempty(stdPred)==1
    stdPred = nan(size(dataCond,1),1);
end

%factors 1,2,3
F.f1.i = dataCond(:,3);
F.f1.nm = 'd';
F.f1.L = unique(F.f1.i);
F.f1.L = sort(F.f1.L,'ascend');
F.f1.n = numel(F.f1.L);

F.f2.i = dataCond(:,2);
F.f2.nm = 'StimStrength';
F.f2.L = unique(F.f2.i);
F.f2.L = sort(F.f2.L,'descend');
F.f2.n = numel(F.f2.L);

F.f3.i = dataCond(:,1);
F.f3.nm = 'Prior std';
F.f3.L = unique(F.f3.i);
F.f3.L = sort(F.f3.L,'descend');
F.f3.n = numel(F.f3.L);

%Graphics
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};

F.f2.colorPre={[0.2 0 0],...
    [0.97 0.2 0],...
    [0.8 0.4 0],...
    [0.3 0.3 0]};

%Mean data and predictions
scrsz=get(0,'ScreenSize');
fig1=figure(1);
set(fig1,'color','w','Position',[0 1.4*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
ymax=nan(F.f2.n,1);
ymin=nan(F.f2.n,1);
h=nan(F.f2.n,1);
count=0;

%graphics
marksz = 100;
ftsz=12;

for j = 1 : F.f2.n
    
    h(j) = subplot(1,F.f2.n,j);
    
    for i = 1 : F.f3.n
        
        %this conditions data
        thisC = F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %data
        count = count+1;
        myPlot1(count) = scatter(F.f1.i( thisC ), meanData( thisC ),marksz,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.f2.color{i},...
            'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));
        
        %pred
        myPlot2 = plot(F.f1.i( thisC ), meanPred( thisC ),...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
    end
    
    %x and ylimits
    ymax(j)=max([meanData;meanPred]);
    ymin(j)=min([meanData;meanPred]);
    xmax(j)=max(F.f1.i);
    xmin(j)=min(F.f1.i);
    
    %x and ylabel
    if j==1
        ylabel('Mean estimates (deg)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)h
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
end

%x and ylimits
set(h,'ylim',[max(ymin) max(ymax)])
set(h,'xlim',[min(xmin) max(xmax)])

%clear up
SLremoveDeadSpace

%selected legend
lg=legend([myPlot1(1:F.f3.n) myPlot2]);
legend(lg,'boxoff')

%-----
%std
%-----
fig2=figure(2);
set(fig2,'color','w','Position',[0 0 scrsz(3)/3 scrsz(4)/3])
ymax=nan(F.f2.n,1);
ymin=nan(F.f2.n,1);
h=nan(F.f2.n,1);
count=0;

for j=1:F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    %axis square
    for i=1:F.f3.n
        
        %this conditions data
        thisC=F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %data
        count=count+1;
        myPlot1(count)=scatter(F.f1.i( thisC ), stdData( thisC ),marksz,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.f2.color{i},...
            'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));
        
        %pred
        myPlot2=plot(F.f1.i( thisC ), stdPred( thisC ),...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linesmoothing','on',...
            'displayName','Bayes');
    end
    
    %x and ylimits
    ymax(j)=max([stdData;stdPred]);
    ymin(j)=min([stdData;stdPred]);
    xmax(j)=max(F.f1.i);
    xmin(j)=min(F.f1.i);
    
    %x and ylabel
    if j==1
        ylabel('Std of estimates (def)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
end

%x and ylimits
set(h,'ylim',[max(ymin) max(ymax)])
set(h,'xlim',[min(xmin) max(xmax)])

%clear up
SLremoveDeadSpace;

%selected legend
lg=legend([myPlot1(1:F.f3.n) myPlot2]);
legend(lg,'boxoff')

%model distribution
function drawDataAndPreDist(Pdata,bins,Ppred,d,StimStrength,pstd,priorModes,...
    cond)

%case meanData or stdData is not input
if isempty(Pdata)==1
    Pdata = nan(size(Ppred));
end
if isempty(Ppred)==1
    Ppred = nan(size(Pdata));
end

%warning
if isempty(Ppred) && isempty(Pdata)==1
    fprintf('%s /n','(drawDataAndPreDist) Both input arguments "Pdata"',...
        ' and "Ppred" are empty. At least one of them should contain data')
end

%we plot the bin value in the middle e.g.., bin 0 to 10 is plotted at
%position 5.
middleBinsForPlot = (bins(2) - bins(1))/2;
bins = bins - middleBinsForPlot;

%graphics
colors = [1 0 0];

%priors, stimStrengths and motion directions as signed linear distance to
%prior mean
Thepriors = unique(pstd);
ThePriorModes = unique(priorModes);
TheStimStrengths = unique(StimStrength);

%motion direction are centered relative to the mean of the experimental
%distributions for both von Mises and bimodal priors.
dlinDisttoPrior = d - nanmean(priorModes,2);

%number of conditions
numPriors = numel(Thepriors);
numStimStrength = numel(TheStimStrengths);

%set of motion directions as signed linear distance to prior mean
%and associated motion direction in deg
[ThedirlinDisttoPrior,pos] = unique(dlinDisttoPrior);
TheDir = d(pos);
motDir = unique(TheDir);
numDir = numel(motDir);

%look at how data distribution change as prior strength increases
%get number of axes and calculate position of plot for each motion
%direction
axesPos = 1:1:numel(ThedirlinDisttoPrior)*numPriors;
dirpos = [];
for i = 1 : numel(Thepriors)
    dirpos = [dirpos; i:numel(Thepriors):numel(axesPos)];
end
dirpos = dirpos';
dirpos = dirpos(:);

%initialize graphics
scrsz = get(0,'ScreenSize');
width = 7;

%initialize the conditions to plot in each figure (StimStrength)
condThisStimStrength = nan(numDir*numPriors,size(cond,2));
condThisStimStrength(:,1) = SLreplicateRows(Thepriors,numDir);
condThisStimStrength(:,3) = repmat(motDir,4,1);
numAllCond = size(condThisStimStrength,1);

%this stimStrength
for j = 1 : numStimStrength
    
    %figure
    figure('color','w',...
        'Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
    
    %this stimStrength
    thisStimStrength = TheStimStrengths(j);
    
    %sync subplots and cond
    condThisStimStrength(:,2) = repmat(TheStimStrengths(j),numDir*numPriors,1);
    for i = 1 : numAllCond
        
        %stop if esc is pressed
        isStop = SLgetKeyPress;
        if isStop == 1 keyboard; end
        
        %time
        t1 = tic;
        
        %axis
        hs(i) = subplot(numDir,numPriors,dirpos(i));
        
        %condition
        thisCon = SLfindRow(condThisStimStrength(i,:),cond);
        thisCon = find(thisCon==1);
        
        %case data
        %draw data and prediction distribution
        if thisCon ~= 0
            hold all
            area(bins,Pdata(:,thisCon),'facecolor',colors,...
                'edgecolor','none')
            plot(bins,Ppred(:,thisCon),'color','k',...
                'linesmoothing','on','linewidth',1);
            
            %graphics
            %dists=[Pdata(:,thisCon); Ppred(:,thisCon)];
            %if sum(isnan(dists))==numel(dists)
            %end
            %ylim is max of all plotted data for ecah axis for clarity.
            %When we use the max of all axis, it's hard to see.
            maxPlot(i) = max([Pdata(:,thisCon); Ppred(:,thisCon)]);
            xlim([0 360])
            ylim([0 maxPlot(i)])
            set(gca,'ytick',[0 fix(maxPlot(i)*10000)/10000],...
                'ytickLabel',[0 fix(maxPlot(i)*10000)/10000])
            
            %motion direction indicator
            plot([cond(thisCon,3) cond(thisCon,3)],[0 0.5*maxPlot(i)],...
                '-','color',[0 0 0],'linewidth',2)
            
            %prior mode
            plot([ThePriorModes ThePriorModes],[0 0.5*maxPlot(i)],...
                '-','color',[0.3 0.6 0.9],'linewidth',2);
            
            %remove axis
            axis off
            set(gca,'ytick',0,'xticklabel',0)
        else
            axis off
        end
        box off
        drawnow
        toc(t1)
        thetime(i)=toc(t1);
        
        %more graphics
        %-------------
        %title top axes
        if sum(dirpos(i)==1:numPriors)==1
            myprior = Thepriors(dirpos(i));
            
            title(['(',num2str(thisStimStrength*100),'% StimStrength -',...
                num2str(myprior),' deg std prior)'],...
                'fontsize',11)
        end
        
        %no x-tick
        set(gca,'xtick',[],'xtickLabel',[],'fontsize',10)
        
        %ylabel left axes
        if sum(dirpos(i)==dirpos(1:numDir))==1
            ylabel('Probability','fontsize',11)
        end
        
        %xlabel bottom axes
        bottomAxes=numAllCond-numPriors+1:numAllCond;
        if sum(dirpos(i)==bottomAxes)==1
            axis on
            xlabel('Estimated directions (deg)','fontsize',11)
            xlim([0 360])
            set(gca,'xtick',25:40:355,'xtickLabel',25:40:355)
        end
    end
    %set(hs,'ylim',[0 max(maxPlot)])
    %clear up
    SLremoveDeadSpace
end

%CMAES
function [fitPtmp,negLogl,counteval,stopflag,outputFit] = fitCMAE(data,...
    feature,...
    stimStrength,...
    pstd,...
    k0,...
    priorShape,...
    priorModes,...
    TheModel,options,varg)

%objective function
%get rid of NaN (cmaes does not like them)
card = k0(8);
if strcmp(TheModel,'withoutCardinal') || isnan(card)
    k0(8)    = [];  %cardinal
    k0(10)   = [];  %tail
    TheModel = 'withoutCardinal';
elseif strcmp(TheModel,'withCardinal') || ~isnan(card)
    k0(11)   = [];  %taik
end

%fit
%note : !!warning !! You code library must have "SLgetLoglBayesianModelcmae.m"
%for cmaes to work
[fitPtmp,negLogl,counteval,stopflag,outputFit] = cmaes('SLgetLoglBayesianModelcmae',k0,sqrt(var(k0)),options,data,feature,stimStrength,pstd,priorShape,priorModes,TheModel,varg);

%rearrange
%when no cardinal parameter
if strcmp(TheModel,'withoutCardinal') || isnan(k0(8))
    %param 11th is nan for now
    %with cmaes
    fitPdum = nan(11,1);
    %likelihood and prior strength
    fitPdum(1:7) = fitPtmp(1:7);
    %cardinal - rand - motor
    fitPdum(9:10) = fitPtmp(8:9);
    fitPtmp = fitPdum;
end
fitPtmp(11) = nan;
fitPtmp = SLmakeRow(fitPtmp);
