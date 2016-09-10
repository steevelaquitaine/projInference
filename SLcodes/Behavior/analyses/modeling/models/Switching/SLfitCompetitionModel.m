

% SLfitCompetitionModel.m
%
%     author: steeve laquitaine
%       date: 131205 (last update 140728)
%    purpose: modeling motion direction estimation behavior with
%             Competition mechanism between learnt priors and evidences.
%
%
%
%
%
%mandatory:
%                                  subjects: e.g., {'sub01'}
%                                    initp : at least 10 params
%
%              'experiment','vonMisesPrior': set von Mises prior data
%                            'bimodalPri15or': or bimodal prior data
%                'vonMisesAndBimodalPriors': combined data
%
%                             'withCardinal': cardinal prior
%                          'withoutCardinal': no cardinal prior
%
%                    'filename','myFilename': set saved file name
%
%analyses
%--------
%
%               - 'MaxLikelihoodFit': fit models to trial data (maximum likelihood)
%                       'fminsearch'
%                            'cmaes'    
%               - 'CrossValR2Fit'   : fit models to data mean and std sorted per
%                                     conditions (cross-validated R-squared)
%                                     (least square fit for Cross-Validation, R^2,
%                                     this requires very good estimates of estimates
%                                     mean and std., i.e., lots of sample per condition)
%
%        - 'CrossValR2FitTrialData': fit models to trial data (least square
%                                    fit for Cross-Validation, R^2)
%
%   - 'CrossValR2FitTrialDataMaxLL': fit models to trial data (max
%                                    likelihood fit for Cross-Validation,logL)
%
%                'modelPredictions': draw data mean and std with models predictions for individual subjects%
%                        'sherlock': run on stanford cluster
%
%
%options
%-------
%
%
%                    'dataPath','/Users/steeve': data path
%                           'filename',filename: saved file
%                                'LoadDataBase': load existing data in dir
%                                       'debug': maxiter = 1
%
%
%
%fitting method
%---------------
%   - 'MaxLikelihoodFit': maximum likelihood fit of trial data
%
%   - 'CrossValR2Fit': least square fit to estimates mean and
%                      std/conditions
%
%   - 'CrossValR2FitTrialData': least square fit to trial data.
%
%   - 'CrossValR2FitTrialDataMaxLL': max likelihood fit to trial data.
%
%
%
%draw data mean and std with models' predictions
%   - "initp" arg: your models parameters
%   - "varargin" :
%           'modelPredictions','withData': show predictions for best fit
%                                          parameters  retrieved
%                                          automatically from path directory
%                                          (if initP are NaN) or input as
%                                          initp.
%           'modelPredictions','noData': draw predictions for input
%                                        parameters without showing data
%'modelPredictions','noData','directoryFitParameter','../modelfit/AIC/','filename','savedModelPred02'
%
%
%-----------
%Description
%-----------
%
%9 fit parameters:
% - 3 von Mises llh concentration parameters k (0<kl<inf, stimulus noise)
% - 4 von Mises (or mixture of von Mises) priors concentration parameters
%   (0<kp<inf)
% - 1 parameter for the fraction of trial with random estimation (0<prand<1)
% - 1 parameter for motor noise (0<km<inf)
%
%
%The model can use:
% - von Mises priors
% - bimodal priors (mixture of von Mises priors)
% -(optionally, cardinal prior, 1 additional fit parameter)
%
%The model can accounts for:
% - cardinal biases (additional Bayesian inference with cardinal priors)
% - random estimation
% - motor noise
%
%references:
%     -Hurliman et al, 2002,VR
%     -Stocker&Simoncelli,2006,NN
%     -Girshick&Simoncelli,2011,NN
%     -Chalk&Series,2012,JoV
%
%
%for examples
%
%       help SLfitCompetitionModelExamples


function [fitP,fitPbkp,R2,sdata,fitPt,negLogl,negLoglbkp,Logl_pertrialBestfit,...
    output] = SLfitCompetitionModel(subjects,...
    initp,...
    varargin)


%call for help
if ieNotDefined('subjects')
    help SLfitCompetitionModel
    return
end

%initialize output
output = [];

%time
t0 = tic;

%backup this .m file in the worspace in output variable
%------------------------------------------------------
%local computer
if ~any(strcmp(varargin,'sherlock'))
    %doesn't work on stanford sherlock cluster
    %     mfilename = SLgetActivemFile;
    %     myMfile = SLbackupMfileInWSpace(mfilename);
    %     output.mfilename  = mfilename;
    %     output.myMfile = myMfile;
end

%set data path
%-------------
%(case von Mises)
%----------------
vrg = varargin;

%local computer
if ~slIsInput(vrg,'sherlock')
    if slIsInput(vrg,'vonMisesPrior')&&~slIsInput(vrg,'bimodalPrior')
        %"dataVM"
        %check and go to input data path
        if slIsInput(vrg,'dataPathVM')
            dataPathVM = varargin{find(strcmp(varargin,'dataPathVM'))+1};
            cd(dataPathVM)
            vararginVM = [varargin,'dataPath',dataPathVM];
            fprintf('%s \n','(SLfitBayesianModel) Creating database ...')
            databankVM = SLMakedatabank(subjects,vararginVM);
            databankBim = [];
            
            %open ui
        elseif slIsInput(vrg,'ui')
            fprintf('%s \n','(SLfitCompetitionModel) Please to set dataPath ...')
            dataPathVM = uigetdir(cd,'Pickup your project e.g., /dataPsychophy/Exp01...');
            cd(dataPathVM)
            vararginVM = [varargin,'dataPath',dataPathVM];
            fprintf('%s \n','(SLfitBayesianModel) Creating database ...')
            databankVM = SLMakedatabank(subjects,vararginVM);
            databankBim = [];
            
            %Load existing database
        elseif slIsInput(vrg,'LoadDataBase')
            if ~isempty(dir('*datbank.mat'))
                fprintf('%s \n','(SLfitBayesianModel) Loading databank in current directory...')
                load('datbank')
                databankVM = databank;
                databankBim = [];
            else
                fprintf('%s \n','(SLfitBayesianModel) Database file datbank.mat does not exist...')
                keyboard
            end
        end
    end
    %cluster
elseif slIsInput(vrg,'sherlock')
    vararginVM  = [varargin,'dataPath','~/data'];

    databankVM  = SLMakedatabank(subjects,vararginVM);
    databankBim = [];
end

%(case bimodal priors)
if slIsInput(varargin,'bimodalPrior') && ~slIsInput(varargin,'vonMisesPrior')
    
    %check and go data path
    if slIsInput(varargin,'dataPathBim')
        dataPathBim = varargin{find(strcmp(varargin,'dataPathBim'))+1};
        cd(dataPathBim)
        
        checkPath = SLexistPath('/Users/steeve_laquitaine/Dropbox/');
        if checkPath == 1
            cd(dataPathBim)
        else
            warning('I couldn t find your path...please enter the path manually.')
            dataPathBim = uigetdir;
        end
    else
        fprintf('%s \n','(SLfitCompetitionModel) Please set dataPathBim')
        return
    end
    databankVM = [];
    vararginBim = [varargin,'dataPath',dataPathBim];
    databankBim = SLMakedatabank(subjects,vararginBim);
end

%(case combined von Mises and bimodal priors)
if slIsInput(varargin,'vonMisesPrior') && ...
        slIsInput(varargin,'bimodalPrior')
    
    %set,go and get databank von Mises prior data path
    if slIsInput(varargin,'dataPathVM')
        dataPathVM = varargin{find(strcmp(varargin,'dataPathVM'))+1};
        cd(dataPathVM)
    else
        fprintf('%s \n','(SLfitCompetitionModel) You need to set dataPathVM')
        return
    end
    vararginVM = [varargin,'dataPath',dataPathVM];
    vararginVM(strcmp(vararginVM,'bimodalPrior'))=[];
    databankVM = SLMakedatabank(subjects,vararginVM);
    
    %set,go and get databank bimodal prior data path.
    if slIsInput(varargin,'dataPathBim')
        dataPathBim = varargin{find(strcmp(varargin,'dataPathBim'))+1};
        cd(dataPathBim)
    else
        fprintf('%s \n','(SLfitCompetitionModel) You need to set dataPathBim')
        return
    end
    vararginBim = [varargin,'dataPath',dataPathBim];
    vararginBim(strcmp(vararginBim,'vonMisesPrior'))=[];
    databankBim = SLMakedatabank(subjects,vararginBim);
end

%for prior modes
%(case von Mises prior)
if slIsInput(varargin,'experiment')
    if slIsInput(varargin,'vonMisesPrior') && ...
            slIsInput(varargin,'bimodalPrior')==0
        priorShape = 'vonMisesPrior';
    end
end

%(case bimodal prior)
if slIsInput(varargin,'experiment')
    if slIsInput(varargin,'bimodalPrior') && ...
            slIsInput(varargin,'vonMisesPrior')==0
        priorShape = 'bimodalPrior';
    end
end

%(case combined von Mises and bimodal prior)
if slIsInput(varargin,'experiment')
    if slIsInput(varargin,'vonMisesPrior') && ...
            slIsInput(varargin,'bimodalPrior')
        priorShape = {'vonMisesPrior','bimodalPrior'};
    end
end

%get the name under which we save the file
if slIsInput(varargin,'filename')
    filename = varargin{find(strcmp(varargin,'filename'))+1};
else
    fprintf('%s \n',['(SLfitCompetitionModel) You need to set the name',...
        ' of the .mat file that will be saved'])
    fprintf('%s \n','(SLfitCompetitionModel) e.g. "filename","mysavedfile"')
end

%(case cardinal prior or not)
if isnan(initp(8))
    TheModel = 'withoutCardinal';
else
    TheModel = 'withCardinal';
end

%Initialize outputs
%(case maximum likelihood fit)
fitP    =[];
R2      =[];
udata   =[];
sdata   =[];
dataOut =[];
predOut =[];
FAbp     =[];
FA       =[];
Sp      =[];
fitPt   =[];
fitPbkp =[];
logL_bkp=[];
logL_pertrialBestfit=[];
negElogL=[];

%(case cross-validated R2)
minSSE=[];
SSE_bkp=[];

%You can set 10 fit parameters (or not,'[]').
fitP = [];

%check that an analysis is input
posTarget = SLfindword(varargin,{'MaxLikelihoodFit',...
    'CrossValR2Fit',...
    'CrossValR2FitTrialData',...
    'modelPredictions'});

if sum(posTarget)==0
    
    %status
    fprintf('%s \n',['(SLfitCompetitionModel) You need to input an analysis:',...
        ' Possible analyses are: '],...
        ' - MaxLikelihoodFit, ',...
        ' - CrossValR2Fit, ',...
        ' - CrossValR2FitTrialData, ',...
        ' - modelPredictions')
end

%(case Maximum likelihood fit)
if slIsInput(varargin,'MaxLikelihoodFit')==1
    
    %status
    fprintf('%s \n','(MLft5) Now fitting the model with maximum likelihood procedure....')
    
    %fit
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,...
        output] = MLft5(databankVM,databankBim,fitP,initp,priorShape,filename,...
        output,varargin);
    
    %save data
    mkdir(['fitComp' subjects{1}])
    cd(['fitComp' subjects{1}])
    slPrintfStr('SLfitCompetitionModel',['Saving data in' pwd])    
    save(filename)
end

%(case Cross-Validated R^2 fit with estimate mean and std per condition)
if slIsInput(varargin,'CrossValR2Fit')==1
    
    %display info
    fprintf('%s \n',['(CrossValR2fit) Now fitting the model with cross-validated R2',...
        ' method....'])
    
    %Coding is quick for now. Ideally I should implement the 3 conditions
    %'vom Mises' only, 'bimodal only' or 'combined' as I did for MLft5
    %function.
    %case von Mises prior only
    if slIsInput(varargin,'vonMisesPrior')
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2fit(databankVM,fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %case bimodal prior only
    if slIsInput(varargin,'bimodalPrior')
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
if slIsInput(varargin,'CrossValR2FitTrialData')==1
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialData) Now fitting the model with',...
        'cross-validated R2 method based on single trial data method....'])
    
    %case von Mises prior only
    if slIsInput(varargin,'vonMisesPrior')
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2FitTrialData(databankVM,...
            fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %case bimodal prior only
    if slIsInput(varargin,'bimodalPrior')
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

%(case Cross-Validated max logLikelihood fit with trial data)
if slIsInput(varargin,'CrossValR2FitTrialDataMaxLL')==1
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Now fitting the model with',...
        'cross-validated max loglikelihood based on single trial data ....'])
    
    %case von Mises prior only
    if slIsInput(varargin,'vonMisesPrior')
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2FitTrialDataMaxLL(databankVM,...
            fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %case bimodal prior only
    if slIsInput(varargin,'bimodalPrior')
        [pred,data,d,coh,pstd,stdPa,output] = CrossValR2FitTrialDataMaxLL(databankBim,fitP,...
            initp,...
            priorShape,...
            filename,...
            output,...
            varargin);
    end
    
    %status
    fprintf('%s \n','(CrossValR2FitTrialDataMaxLL) Ending Cross-validation.')
    
    %save data
    save(filename)
    fprintf('%s \n','(CrossValR2FitTrialDataMaxLL) .mat file is saved.')
end

% %(Case draw data mean, std and distribution with their best models'
% %predictions)
% %----------------------------------------------------------------
% if slIsInput(varargin,'modelPredictions')==1
%
%     %case von Mises prior only
%     if slIsInput(varargin,'vonMisesPrior') && ...
%             slIsInput(varargin,'bimodalPrior')==0
%         databank = databankVM;
%     end
%
%     %case bimodal prior only
%     if slIsInput(varargin,'bimodalPrior') &&...
%             slIsInput(varargin,'vonMisesPrior')
%         databank = databankBim;
%     end
%
%     %case we show data
%     %-----------------
%     if slIsInput(varargin,'withData'))==1
%
%         %display info
%         fprintf('%s \n',['(drawMeanPreCentered) Now draw data mean,'...
%             'std and distributions with their best models predictions'])
%
%
%         %number of subjects
%         subjects = unique(databank.subjects);
%         numSub = numel(subjects);
%
%         %warn that data and predictions will be averaged across subjects
%         if numSub > 1
%             fprintf('%s \n',['(drawMeanPreCentered) More that 1 subject have been found',...
%                 '. Data and models predictions will be averaged across subjects'])
%         end
%
%         %all experimental conditions
%         output.uniqCond = SLuniqpair([databank.Pstd databank.coherence ...
%             databank.motiondirDeg]);
%         numcond = size(output.uniqCond,1);
%
%         %preallocate
%         output.meanData = nan(numcond,numSub);
%         output.stdData  = nan(numcond,numSub);
%         output.meanPred = nan(numcond,numSub);
%         output.stdPred  = nan(numcond,numSub);
%
%         %get subject's best fit parameters
%         %data come from modelfit folder. The tree structure of the folders
%         %within the  experiment folder must remain the same. The name of
%         %the folder must not change
%         %!! can probably improve here by directly getting all data from
%         %modelfit folder.
%         %cd ../modelfit/model_CompDiv/
%         dirModelFit = dir('../modelfit/AIC/model_CompDiv/sub*');
%         subjectsInDir = {dirModelFit.name};
%         numSubInDir = numel(subjectsInDir);
%         subjectsInDirlist = nan(numSubInDir,1);
%
%         for ij = 1 : numel(subjectsInDir)
%             subjectsInDirlist(ij) = str2double(subjectsInDir{ij}(4:end));
%         end
%
%         %data and predictions for each subject
%         for sub = 1 : numSub
%
%             %status
%             fprintf('%s \n',['(drawMeanPre) Retrieving data for subject ',...
%                 num2str(sub)])
%
%             %this trial's subject
%             TrialsThisSub = subjects(sub)==databank.subjects;
%
%             %get data and conditions
%             data = round(databank.estimatesDeg(TrialsThisSub));
%             d = databank.motiondirDeg(TrialsThisSub);
%             coh = databank.coherence(TrialsThisSub);
%             pstd = databank.Pstd(TrialsThisSub);
%             priorModes = cell2mat(databank.priormodes(TrialsThisSub,:));
%
%             %data mean and std
%             [meanDatatmp,stdDatatmp,dataCondtmp] = makeDataMeanAndStd(data,...
%                 d,coh,pstd,priorModes,priorShape);
%
%             %Case we want best fit parameters, show data and predictions
%             %-----------------------------------------------------------
%             %Best fit parameters from maximum likelihood fit are automatically
%             %loaded from the .../modelfit/AIC folder.
%             %they are automatically loaded
%             if sum(isnan(initp))==length(initp)
%
%                 %find subjectsInDir that matches subjects append this subject to
%                 %loading directory, find datafit mat file, append to directory and
%                 %load model's best fit parameters
%                 thisSubInDir = subjectsInDir(subjectsInDirlist==subjects(sub));
%                 FitfilenameThisSub = dir(['../modelfit/AIC/model_CompDiv/' ...
%                     cell2mat(thisSubInDir) '/datafit*']);
%                 FitfilenameThisSub = FitfilenameThisSub.name;
%                 load(['../modelfit/AIC/model_CompDiv/' thisSubInDir{1},'/',...
%                     FitfilenameThisSub(1:end-4)],'fitP')
%
%                 %store subjects' best fit parameters
%                 output.fitP(sub,:) = fitP.p;
%                 clear fitP
%             end
%
%             %case simulation parameters are input (free simulation)
%             %------------------------------------------------------
%             if sum(isnan(initp))~=length(initp)
%                 output.fitP(sub,:) = initp;
%             end
%
%             %models' estimate mean, std and distribution predictions
%             [meanPredtmp,stdPredtmp,cond,PestimateGivenModelUniq,...
%                 MAP,output] = SLmakePredictionsCompetitionModel(d,coh,pstd,output.fitP(sub,:),priorShape,...
%                 priorModes,[],output,varargin);
%
%             %subject's estimate distribution
%             %[pData,xpdf] = makeDataDist(data,d,coh,pstd,priorModes);
%             [pData,xpdf] = makeDataDist(data,d,coh,pstd,priorModes,cond);
%
%             %make sure predicted and data distributions are calculated on
%             %the same space. Works for predicted estimate space = 1:1:360 deg.
%             %adjust predicted estimate probability distribution from 1:1:360 to
%             %0:10:360 by summing probabilities within consecutive bins of 10
%             %deg (law of probabilities).
%             commonSpace = xpdf;
%             numSpace = numel(commonSpace)-1;
%             [~,bins]= histc(0:1:360,commonSpace);
%             bins(end) = [];
%             bins = bins';
%             numCond = size(cond,1);
%             pPred = nan(numSpace,numCond);
%             for ijk = 1 : numCond
%                 pPred(:,ijk) = SLcumSumInBin(PestimateGivenModelUniq(:,ijk),bins);
%             end
%             commonSpace = commonSpace(2:end);
%
%             %match conditions, data and predictions across predictions and
%             %data and across subjects
%             numPredEst = size(pPred,1);
%             numDataEst = size(pData,1);
%             for j = 1 : size(output.uniqCond,1)
%
%                 %case cond (predictions) and dataCondtmp (data) are same as
%                 %expected
%                 %get current condition sync with a reference ordering of the
%                 %conditions
%                 if isequal(dataCondtmp,cond)
%                     [~,posThisCond] = ismember(output.uniqCond(j,:),...
%                         dataCondtmp,'rows');
%
%                     %case condition exists
%                     if ~isequal(posThisCond,0)
%
%                         %match data and predictions mean, std and distributions
%                         %between conditions and subjects to the order in
%                         %output.uniqCond
%                         output.meanData(j,sub) = meanDatatmp(posThisCond);
%                         output.stdData(j,sub)  = stdDatatmp(posThisCond);
%                         output.meanPred(j,sub) = meanPredtmp(posThisCond);
%                         output.stdPred(j,sub)  = stdPredtmp(posThisCond);
%                         output.PdisPred(:,j,sub) = pPred(:,posThisCond);
%                         output.PdisData(:,j,sub) = pData(:,posThisCond);
%                     else
%                         %else NaN
%                         output.meanData(j,sub) = NaN;
%                         output.stdData(j,sub)  = NaN;
%                         output.meanPred(j,sub) = NaN;
%                         output.stdPred(j,sub)  = NaN;
%                         output.PdisPred(:,j,sub) = NaN(numPredEst,1);
%                         output.PdisData(:,j,sub) = NaN(numDataEst,1);
%                     end
%                 end
%             end
%         end
%
%         %average over subjects
%         %---------------------
%         %data mean, std and distributions
%         for i = 1 : size(output.meanData,1)
%             CircMeanOvSub = SLcircMeanStd(output.meanData(i,:)','polar');
%             output.meanDataOvSub(i) = CircMeanOvSub.deg.mean;
%         end
%         output.meanDataOvSub = SLmakeColumn(output.meanDataOvSub);
%         output.stdDataOvSub     = nanmean(output.stdData,2);
%         output.meanDisDataOvSub = nanmean(output.PdisData,3);
%
%         %predictions mean and std and distributions
%         for i = 1 : size(output.meanPred,1)
%             PredCircMeanOvSub = SLcircMeanStd(output.meanPred(i,:)','polar');
%             output.meanPredOvSub(i) = PredCircMeanOvSub .deg.mean;
%         end
%         output.meanPredOvSub = SLmakeColumn(output.meanPredOvSub);
%         output.stdPredOvSub     = nanmean(output.stdPred,2);
%         output.meanDisPredOvSub = nanmean(output.PdisPred,3);
%
%         %draw data mean and std and their models' predictions
%         drawMeanPreCentered(output.meanDataOvSub,output.stdDataOvSub,...
%             output.meanPredOvSub,output.stdPredOvSub,...
%             output.uniqCond)
%
%         %Code predictions about estimate distribution here
% %         drawDataAndPreDistCentered(output.meanDisDataOvSub,commonSpace,...
% %             output.meanDisPredOvSub,...
% %             d,coh,pstd,priorModes,output.uniqCond,varargin)
%
%         SLdrawModelsPredictionHistwithDataCentered(output.meanDisDataOvSub,commonSpace,...
%             output.meanDisPredOvSub,d,coh,pstd,priorModes,output.uniqCond,...
%             priorShape,varargin)
%     end
%
%     %case show data
%     %--------------
%     if slIsInput(varargin,'noData')==1
%         %Case parameters are input and we want to simulate predictions without
%         %showing data. Works with many subjects too.
%         %---------------------------------------------------------------------
%         if sum(isnan(initp))~=length(initp)
%
%             %status
%             fprintf('%s \n','(drawMeanPreCentered) Simulating model predictions...')
%
%             %all experimental conditions
%             d = databank.motiondirDeg;
%             pstd = databank.Pstd;
%             coh = databank.coherence;
%             priorModes = cell2mat(databank.priormodes);
%             output.uniqCond = SLuniqpair([pstd coh d]);
%             numcond = size(output.uniqCond,1);
%
%             %models' estimate mean, std and distribution predictions
%             output.fitP = initp;
%             [meanPred,stdPred,cond,PestimateGivenModelUniq,...
%                 MAP] = SLmakePredictionsCompetitionModel(d,coh,pstd,output.fitP,priorShape,...
%                 priorModes,[],output,varargin);
%
%             %make sure predicted and data distributions are calculated on
%             %the same space. Works for predicted estimate space = 1:1:360 deg.
%             %adjust predicted estimate probability distribution from 1:1:360 to
%             %0:10:360 by summing probabilities within consecutive bins of 10
%             %deg (law of probabilities).
%             %commonSpace = 0:10:360;
%             commonSpace = 0:10:360;
%             numSpace = numel(commonSpace)-1;
%             [~,bins] = histc(0:1:360,commonSpace);
%             bins(end) = [];
%             bins = bins';
%             numCond = size(cond,1);
%             pPred = nan(numSpace,numCond);
%             for ijk = 1 : numCond
%                 pPred(:,ijk) = SLcumSumInBin(PestimateGivenModelUniq(:,ijk),bins);
%             end
%             commonSpace = commonSpace(2:end);
%
%             %predicted mean,std and distributions
%             output.meanPred    = meanPred;
%             output.stdPred     = stdPred;
%             output.meanDisPred = pPred;
%             output.uniqCond    = cond;
%
%             %draw predicted mean, std and distribution
%             %mean and std
%             drawMeanPreCentered([],[],output.meanPred,output.stdPred,output.uniqCond)
%
%             %distribution
% %             drawDataAndPreDistCentered([],commonSpace,...
% %                 output.meanDisPred,...
% %                 d,coh,pstd,priorModes,output.uniqCond,varargin)
%
%             SLdrawModelsPredictionHistwithDataCentered([],commonSpace,...
%                 output.meanDisPred,...
%                 d,coh,pstd,priorModes,output.uniqCond,priorShape,varargin)
%         end
%     end
% end

%(Case data mean, std and distribution with models'predictions)
if slIsInput(varargin,'modelPredictions')==1
    
    %case von Mises prior
    if slIsInput(varargin,'vonMisesPrior') && slIsInput(varargin,'bimodalPrior')==0
        databank = databankVM;
    end
    
    %case bimodal prior
    if slIsInput(varargin,'bimodalPrior') && slIsInput(varargin,'vonMisesPrior')==0
        databank = databankBim;
    end
    
    %plot data or not
    if slIsInput(varargin,'withData')==0 && slIsInput(varargin,'noData')==0
        
        %please input with or without data
        error('(SLfitCompetitionModel) Please input "withData" or "noData" argument to display or not the data...')
        keyboard
        
    end
    
    %with data
    if slIsInput(varargin,'withData')==1
        
        %status
        fprintf('%s \n',['(SLdrawModelsPredictionCentered) Now draw data mean,'...
            'std and distributions with best models predictions'])
        
        %model(cardinal or not)
        %without cardinal
        if isnan(initp(8))
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Your model has no cardinal prior.')
            modelFold='model_CompDiv';
            
            %with cardinal
        elseif ~isnan(initp(8))
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Your model has cardinal priors')
            modelFold = 'model_CompDiv_withCard';
        end
        
        %number of subjects
        inputSubjects = subjects;
        subjects = unique(databank.subjects);
        numSub = numel(subjects);
        
        %warn that data and predictions will be averaged across subjects
        if numSub > 1
            fprintf('%s \n',['(SLdrawModelsPredictionCentered) More that 1 subject have been found',...
                '. Data and models predictions will be averaged across subjects'])
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
                output.uniqCond = SLuniqpair([priorCond databank.stimStrength ...
                    databank.stimFeatureDeg]);
                numcond = size(output.uniqCond,1);
            end
            
            %Warning if priorShape is missing
        else
            fprintf('%s \n',['(SLfitCompetitionModel : ModelPredictions)',...
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
        if slIsInput(varargin,'directoryFitParameter')
            
            %directory
            inputdir = varargin{find(strcmp(varargin,...
                'directoryFitParameter')) + 1};
            dirModelFit = dir([inputdir,modelFold,'/sub*']);
            
            %Check that the directory exists
            if isempty(dirModelFit)==1
                fprintf('%s \n',['(SLfitCompetitionModel) !! WARNING !! Sorry but',inputdir,modelFold,'/sub does not exist.....'])
                keyboard
            end
            
            %subjects in directory
            subjectsInDir = {dirModelFit.name};
            numSubInDir = numel(subjectsInDir);
            
            %status
            fprintf('%s \n','(SLfitCompetitionModel) Directory for best fit',...
                '  parameters was input manually')
            
            %case we use input Parameters
            %----------------------------
        elseif slIsInput(varargin,'inputFitParameters')
            inputdir = initp;
            numSubInDir = numSub;
            subjectsInDir = inputSubjects;
            
            %case not input
            %--------------
        else
            fprintf('%s \n','(SLfitCompetitionModel) Please input best fit params directory after "directoryFitParameter"...')
            fprintf('%s \n','(SLfitCompetitionModel) e.g.,"directoryFitParameter","../modelfit/AIC/"')
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
            fprintf('%s \n',['(SLdrawModelsPredictionCentered) Retrieving data for subject ',...
                cell2mat(subjectsInDir)])
            
            %this trial's subject
            TrialsThisSub = subjects(sub)==databank.subjects;
            
            %data and conditions
            data = round(databank.estimatesDeg(TrialsThisSub));
            d    = databank.stimFeatureDeg(TrialsThisSub);
            coh  = databank.stimStrength(TrialsThisSub);
            pstd = databank.Pstd(TrialsThisSub);
            priorModes = cell2mat(databank.priormodes(TrialsThisSub,:));
            
            %data mean and std
            [meanDatatmp,stdDatatmp,dataCondtmp] = makeDataMeanAndStd(data,...
                d,coh,pstd,priorModes,priorShape);
            
            %Case we want best fit parameters, show data and predictions
            %-----------------------------------------------------------
            %Directory for best fit parameters from maximum likelihood fit
            %must be input (e.g.,.../modelfit/..AIC folder).
            if sum(isnan(initp))==length(initp)
                
                %find subjectsInDir that matches subjects append this subject to
                %loading directory, find datafit mat file, append to directory and
                %load model's best fit parameters
                if slIsInput(varargin,'directoryFitParameter')
                    
                    %subject
                    thisSubInDir = subjectsInDir(subjectsInDirlist==subjects(sub));
                    
                    %file
                    FitfilenameThisSub = dir([inputdir,modelFold,'/' ...
                        cell2mat(thisSubInDir) '/datafit*']);
                    FitfilenameThisSub = FitfilenameThisSub.name;
                    
                    %load
                    load([inputdir,modelFold,'/',thisSubInDir{1},'/',...
                        FitfilenameThisSub(1:end-4)],'fitP')
                    
                    %Warning
                else
                    fprintf('%s \n',['(SLfitCompetitionModel) You need to input',...
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
                clear fitP
            end
            
            %case simulation parameters are input (free simulation)
            %------------------------------------------------------
            if sum(isnan(initp))~=length(initp)
                output.fitP(sub,:) = initp;
                
                %status
                fprintf('%s \n','(SLfitCompetitionModel) I am using the input model simulation parameters')
            end
            
            %models' estimate mean, std and distribution predictions
            [meanPredtmp,stdPredtmp,cond,PestimateGivenModelUniq,...
                MAP] = SLmakePredictionsCompetitionModel(d,coh,pstd,output.fitP(sub,:),priorShape,...
                priorModes,[],output,varargin);
            
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
        
        %Code predictions about estimates' distributions here
        SLdrawModelsPredictionHistwithDataCentered(output.meanDisDataOvSub,commonSpace,...
            output.meanDisPredOvSub,d,coh,pstd,priorModes,output.uniqCond,...
            priorShape,varargin)
    end
    
    %no data
    if slIsInput(varargin,'noData')==1
        
        %Case input param & predictions without data
        if sum(isnan(initp))~=length(initp)
            
            %status
            fprintf('%s \n','(SLdrawModelsPredictionCentered) Simulating model predictions...')
            
            %status: model
            if ~isnan(initp(8))
                fprintf('%s \n','(SLdrawModelsPredictionCentered) Model has cardinal prior')
            end
            if isnan(initp(8))
                fprintf('%s \n','(SLdrawModelsPredictionCentered) Model has no cardinal prior')
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
                MAP] = SLmakePredictionsCompetitionModel(d,...
                coh,...
                pstd,...
                output.fitP,...
                priorShape,...
                priorModes,...
                [],...
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
            
            %distribution
            SLdrawModelsPredictionHistwithDataCentered([],commonSpace,...
                output.meanDisPred,...
                d,coh,pstd,priorModes,output.uniqCond,priorShape,varargin)
        end
    end
end

%Case get standard error for best fit parameters
if slIsInput(varargin,'stdBestfitP')==1
    
    %status
    fprintf('%s \n',['(MLft5) Now fitting the model with',...
        'maximum likelihood method to get best fit parameters standard error'])
    
    %fit
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,...
        output] = MLft5(databank,fitP,initp,priorShape,...
        filename,output,varargin);
    
    %save data
    save(filename)
end

%time
output.fitDuration = toc(t0);














%factors
function FA = initFAs(databank,FAs)

%status
fprintf('%s \n','(initFAs) Now getting factors...')

%data
est_data=cell2mat(databank.data(:,(strcmp(databank.nm,'est_coor'))==1));
es.dir=cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1));

%Set factors
%FA 1
FA.g1.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{1}  ))==1));
FA.g1.nm       =FAs{1};
%FA 2
FA.g2.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{2}   ))==1));
FA.g2.nm       =FAs{2};
%FA 3
FA.g3.thisT=cell2mat(databank.data(:,(strcmp(databank.nm,FAs{3}  ))==1));
FA.g3.nm       =FAs{3};


%CASE 1: Prior is unimodal
%-------------------------
%Get and order levels of each group
%FA1
FA.g1.lvlsnm=unique(FA.g1.thisT);
FA.g1.lvlsnm=sort(FA.g1.lvlsnm,'descend');
FA.g1.lvlsnb=numel(FA.g1.lvlsnm);
clear i
for i=1:FA.g1.lvlsnb
    index.g1.lvli(i)={find(FA.g1.thisT==FA.g1.lvlsnm(i))};
end

%FA2
FA.g2.lvlsnm=unique(FA.g2.thisT);
FA.g2.lvlsnm=sort(FA.g2.lvlsnm,'descend');
FA.g2.lvlsnb=numel(FA.g2.lvlsnm);
for i=1:FA.g2.lvlsnb
    index.g2.lvli(i)={find(FA.g2.thisT==FA.g2.lvlsnm(i))};
end

%FA3
FA.g3.lvlsnm=unique(FA.g3.thisT);
FA.g3.lvlsnm=sort(FA.g3.lvlsnm,'ascend');
FA.g3.lvlsnb=numel(FA.g3.lvlsnm);
for i=1:FA.g3.lvlsnb
    index.g3.lvli(i)={find(FA.g3.thisT==FA.g3.lvlsnm(i))};
end

%-----------------------
%Case max likelihood fit
%-----------------------
function [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
    MLft5(databankVM,databankBim,fitPinput,...
    initp,...
    priorShape,...
    filename,...
    output,...
    varargin)

%status
fprintf('%s \n','(SLfitCompetitionModel) Fitting trial-data (maximum likelihood)...')

%(case check fitting fake data).
%Get simulated data present in the workspace.
%-------------------------------------------
if sum(strcmp(varargin{1},'testFakedata'))
    
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
    
    %databank
    if ~isempty(databankVM)
        databank = databankVM;
    elseif isempty(databankBim)
        databank = databankBim;
    end
    
    %subjects' data
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'estimatedFeature'))==1)));
    
    %task conditions
    disp         = cell2mat(databank.data(:,(strcmp(databank.nm,'FeatureSample'))==1));
    stimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd         = cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes   = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make 360 and 0 degrees same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    stimStrength  = stimStrength(isnan(data)==0);
    pstd = pstd(isnan(data)==0);
    priorModes = priorModes(isnan(data)==0,:);
    
    %store
    output.data = data;
    output.disp = disp;
    output.stimStrength = stimStrength;
    output.pstd = pstd;
    output.priorModes = priorModes;
end

%(case combined von Mises and bimodal prior)
%-------------------------------------------
if ~isempty(databankVM) && ~isempty(databankBim)
    
    %status
    fprintf('%s \n','(MLft5) Fitting data for combined von Mises and bimodal prior experiments...')
    
    %(case von Mises prior)
    %----------------------
    %status
    fprintf('%s \n','(MLft5) Collecting von Mises prior data')
    
    %subjects' data
    dataVM = round(cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'est_dir'))==1)));
    
    %experimental conditions
    dispVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'FeatureSample'))==1));
    stimStrengthVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'stimStrength'))==1));
    pstdVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'Pstd'))==1));
    priorModesVM = cell2mat(databankVM.data(:,(strcmp(databankVM.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    dataVM(dataVM==0) = 360;
    
    %remove missing data
    dataVM = dataVM(isnan(dataVM)==0);
    dispVM = dispVM(isnan(dataVM)==0);
    stimStrengthVM  = stimStrengthVM(isnan(dataVM)==0);
    pstdVM = pstdVM(isnan(dataVM)==0);
    priorModesVM = priorModesVM(isnan(dataVM)==0,:);
    
    %store
    output.dataVM = dataVM;
    output.dispVM = dispVM;
    output.stimStrengthVM = stimStrengthVM;
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
    stimStrengthBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'stimStrength'))==1));
    pstdBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'Pstd'))==1));
    priorModesBim = cell2mat(databankBim.data(:,(strcmp(databankBim.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    dataBim(dataBim==0) = 360;
    
    %remove missing data
    dataBim = dataBim(isnan(dataBim)==0);
    dispBim = dispBim(isnan(dataBim)==0);
    stimStrengthBim  = stimStrengthBim(isnan(dataBim)==0);
    pstdBim = pstdBim(isnan(dataBim)==0);
    priorModesBim = priorModesBim(isnan(dataBim)==0,:);
    
    %store
    output.dataBim = dataBim;
    output.dispBim = dispBim;
    output.stimStrengthBim = stimStrengthBim;
    output.pstdBim = pstdBim;
    output.priorModesBim = priorModesBim;
end

%for fminsearch - Nelder-Mead
%OptionsNM=optimset('Display','iter');
OptionsNM=[];

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
    kllh       = initp(1:3);
    klearnt    = initp(4:7);
    kcardinal  = initp(8);
    fractRand  = initp(9);
    motorNoise = initp(10);
    
    %(Case we want to fit data with maximum likelihood, fminsearch)
    %-------------------------------------------------------------
    if sum(strcmp(varargin{1},'MaxLikelihoodFit'))==1
        
        %(1) Input parameters (best matching parameters)
        k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise];
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
        k0(10,:) = [repmat(nanmean(k0(1,1:8)),1,8) fractRand motorNoise];
    end
    
    %(Case we want to get best fit parameters' standard error, fmincon)
    %------------------------------------------------------------------
    %Only keep input parameters (fit parameters e.g., from a previous fit
    %for examples)
    if sum(strcmp(varargin{1},'stdBestfitP'))==1
        
        %status
        fprintf('%s \n',['(MLft5) Setting your input initial parameters as'...
            'the only set of initial parameter'])
        k0(1,:) = [kllh klearnt kcardinal fractRand motorNoise];
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
    
    %(Case we fit data with max likelihood
    %-------------------------------------
    if sum(strcmp(varargin{1},'MaxLikelihoodFit'))==1
        
        %(case von Mises or Bimodal prior)
        %---------------------------------
        if isempty(databankVM) || isempty(databankBim)
            
            %fitting options
            %---------------
            %debugging
            %options = optimset('MaxIter',1,'MaxFunEvals',1);
            options = optimset('MaxIter',200*size(k0,1),...
                'MaxFunEvals',10*200*size(k0,1));
            %options = [];
            
            if any(strcmp(varargin{:},'fminsearch'))
                parfor i = 1 : size(k0,1)
                    %for i = 1 : size(k0,1)
                    
                    %status
                    fprintf('%s \n',['(MLft5)',' Set of initial parmeters: ',...
                        num2str(i),'/',num2str(size(k0,1))])
                    
                    %Nelder-Mead
                    [fitPtmp,negLogl,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                        getLogl(data,...
                        disp,...
                        stimStrength,...
                        pstd,...
                        fitPtmp,...
                        priorShape,...
                        priorModes,...
                        TheModel,varargin{1}),...
                        k0(i,:),...
                        options);
                    
                    %Fit parameters and SSE
                    fitPbkp(i,:) = fitPtmp;
                    negLoglbkp(i) = negLogl;
                    exitflagbkp{i} =  exitflag;
                    outputFitbkp{i} =  outputFit;
                end
                output.fitAlgo = 'fminsearch';
                
                %covariance matrix a
            elseif any(strcmp(varargin{:},'cmaes'))
                
                %options
                options = cmaes;
                options.LBounds = 0;  %no < 0
                options.TolFun = 1e-4; %stop when obj(t2) - obj(t1) < 0.5
                
                %debug mode
                if slIsInput(varargin{:},'debug')
                    fprintf('%s \n','(SLfitCompetitionModel) Debug: MaxIter = 1.')
                    options.MaxIter = 1;
                end

                %testing with defined maxiter for speed----
                options.MaxIter = 200*9;
                options.MaxFunEvals = options.MaxIter*10;
                %------------------------------------------
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
                    fprintf('\n',['(MaxLikelihoodFit)',' setting initial parameters for cmaes: ',num2str(i),'/',num2str(size(k0,1))])
                    
                    %fit
                    [fitPtmp,negLogl,~,stopflag,outputFit] = fitCMAE(data,disp,stimStrength,pstd,k0(i,:),priorShape,priorModes,TheModel,options,varargin{1});
                    
                    %Fit parameters and SSE
                    fitPbkp(i,:) = fitPtmp;
                    negLoglbkp(i) = negLogl;
                    exitflagbkp{i} = stopflag;
                    outputFitbkp{i} = outputFit;
                end
                output.fitAlgo = 'cmaes';
                output.fitoptions = options;
            end
            
            %store
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
                [~,~,Logl_pertrialBestfit] = getLogl(data,disp,stimStrength,pstd,fitP.p,...
                    priorShape,priorModes,TheModel,varargin{1});
                
                %cmaes fitting
            elseif any(strcmp(varargin{:},'cmaes'))
                
                %no cardinal
                if strcmp(TheModel,'withoutCardinal') || isnan(fitP.p(8))
                    fitP.p(8)  = [];
                    TheModel = 'withoutCardinal';
                end
                [~,~,Logl_pertrialBestfit] = SLgetLoglCompCmae(fitP.p',data,disp,stimStrength,pstd,...
                    priorShape,priorModes,TheModel,varargin{1});
                %rearrange
                if strcmp(TheModel,'withoutCardinal') || isnan(fitP.p(8))
                    fitPdum       = nan(10,1); %card
                    fitPdum(1:7)  = fitP.p(1:7);
                    fitPdum(9:10) = fitP.p(8:9);
                end
                fitP.p = fitPdum;
            end
            %store
            output.Logl_pertrialBestfit = Logl_pertrialBestfit;
            output.fitP = fitP.p;
        end
        
        %(case combined von Mises and Bimodal prior)
        %-------------------------------------------
        if ~isempty(databankVM) && ~isempty(databankBim)
            
            %fitting options
            %---------------
            %(debugging)
            %options = optimset('MaxIter',1,'MaxFunEvals',1);
            options = [];
            
            %fit
            %parfor i = 1 : size(k0,1)
            for i = 1 : size(k0,1)
                
                %status
                fprintf('%s \n',['(MLft5)',' Set of initial parmeters: ',...
                    num2str(i),'/',num2str(size(k0,1))])
                
                %Nelder-Mead
                [fitPtmp,negLogl,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                    getLoglCombinedPriors(dataVM,dataBim,...
                    dispVM,dispBim,...
                    stimStrengthVM ,stimStrengthBim,...
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
                stimStrengthVM ,stimStrengthBim,...
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
        
        %parfor i = 1 : size(k0,1)
        for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n','(MLft5) Set of the initial parameters...')
            
            %fmincon
            [fitPtmp,negLogl,exitflag,outputFit,LAMBDA,...
                GRAD,fitHessian] = fmincon( @(fitPtmp) ...
                getLogl(data,...
                disp,...
                stimStrength,...
                pstd,...
                fitPtmp,...
                priorShape,...
                priorModes,...
                TheModel,varargin{1}),...
                k0(i,:),...
                [],[],[],[],...
                [0       0    0    0    0    0    0                0 0    0],...
                [3000 3000 3000 3000 3000 3000 3000 higherConstKcard 1 3000],...
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
        
        [~,~,Logl_pertrialBestfit] = getLogl(data,disp,stimStrength,pstd,fitP.p,...
            priorShape,priorModes,TheModel,varargin{1});
        
        %store
        output.Logl_pertrialBestfit = Logl_pertrialBestfit;
        output.fitP = fitP.p;
    end
    
    
    %(case von Mise prior)
    %---------------------
    if sum(strcmp(priorShape,'vonMisesPrior')) && ...
            ~sum(strcmp(priorShape,'bimodalPrior'))
        fitP.nm = {'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10',...
            'kcardinal','Prand','km'};
    end
    
    %(case bimodal prior)
    %---------------------
    if sum(strcmp(priorShape,'bimodalPrior')) && ...
            ~sum(strcmp(priorShape,'vonMisesPrior'))
        fitP.nm = {'coh24','coh12','coh6','pstd145_305','pstd165_285',...
            'pstd185_265','pstd205_245',...
            'kcardinal','Prand','km'};
    end
    
    %(case von Mises and bimodal prior)
    %----------------------------------
    if sum(strcmp(priorShape,'vonMisesPrior')) && ...
            sum(strcmp(priorShape,'bimodalPrior'))
        fitP.nm = {'coh24','coh12','coh6','pstd80modes225or20modes145_305',...
            'pstd40modes225or20modes165_285',...
            'pstd20modes225or20modes185_265',...
            'pstd10modes225or20modes205_245',...
            'kcardinal','Prand','km'};
    end
    
    %(case no cardinal prior)
    %------------------------
    %Kcardinal=NaN.
    if isnan(kcardinal)
        fitP.p(strcmp(fitP.nm,'kcardinal')) = NaN;
    end
    
    %parameters std
    if sum(strcmp(varargin{1},'stdBestfitP'))==1
        [~,~,~,~,~,output] = SLmakePredictionsCompetitionModel(disp,stimStrength,pstd,fitP.p,priorShape,...
            priorModes,'Trial',output,varargin{1});
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
    
    [negLogl,~,Logl_pertrialBestfit] = getLogl(data,disp,stimStrength,pstd,fitPinput,...
        priorShape,priorModes,TheModel,varargin{1});
    negLoglbkp = negLogl;
    fitPbkp = fitPinput;
    fitP.p = fitPinput;
    fitP.nm = {'coh24','coh12','coh6','klearnt80','klearnt40','klearnt20',...
        'klearnt10','Kcardinal','Prand','km'};
    
    %store data
    save(filename)
end

%logl
function [negLogl,fitP,Logl_pertrial] = getLogl(data,d,stimStrength,pstd,fitP,...
    priorShape,priorModes,TheModel,varargin)
%Make predictions
%----------------
%   Data from each trial can be explained by either:
%
%       Bayesian inference with sensory evidence and cardinal
%       -----------------------------------------------------
%       measurement mi (=llh mean: ul) is drawn from a von Mises measurement
%       density (v(mi,km)) in 1-("Pp"+"Prandom") fraction of trials.
%       is converted to likelihood and combined with learnt prior. The
%       probability of all possible percepts is the probability that it
%       comes from a given measuremen mi given by (v(mi,km)) * the
%       probability that subject switch given by a divisive competition
%       mechanism between evidence strength and learnt prior strength.
%
%       Bayesian inference with learnt prior and cardinal prior
%       ---------------------------------------------------------
%       The probability of a percept given subject switched to the learnt
%       prior is the probability to switch given by the divisive
%       competition * the probability that this percept is produced by the
%       model(e.g., if cardinal is flat P(learnt prior mean)=100%;  or if
%       learnt prior is flat P(each of 4 cardinal modes) = 25%)
%
%       Random estimation
%       -----------------
%       an estimate is produced randomly in a fraction of trials ("Prandom").
%
%       And motor noise
%       ---------------
%       ~v(0,kmo) is added to the estimate (let's call it percept)
%       of those two processes.
%
%       Thus, the maximum likelihood of observing each trial data assuming this
%       model is:
%       probability of data given Bayesian inference * (1-"Pp"+"Prandom") +
%       probability of data given random choice * "Prandom".
%       then both are convolved with motor noise (note: its is because convolution
%       is distributive).

%get fit parameters. They should be in this order.
%(von Mises prior)
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kcardinal','Prand','km'
%(bimodal prior)
%'coh24','coh12','coh6','[145 305]','[165 285]','[185 265]','[205 245]','kcardinal','Prand','km'
%time
ticLogL = tic;

%model parameters
kl1 = fitP(1);
kl2 = fitP(2);
kl3 = fitP(3);
klearnt1 = fitP(4);
klearnt2 = fitP(5);
klearnt3 = fitP(6);
klearnt4 = fitP(7);
kcardinal = fitP(8);
Prandom = fitP(9);
km = fitP(10);
ulall = 1:1:360;

%--------------
%SET THE PRIORS
%--------------
%
%(case cardinal prior)
%---------------------
%case we don't want to fit the cardinal prior.
if isnan(kcardinal)
    TheModel='withoutCardinal';
    kcardinal=0;
else
    %case we don't want to fit the cardinal prior.
    TheModel='withCardinal';
end

%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time.
%Constrain other fraction parameters between 0 and 1, and other parameters
%as>0.
%and case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom > 1
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    return
end
if any(fitP < 0)
    negLogl = 1e09;
    Logl_pertrial= - 1e09;
    return
end
if any(isnan(fitP([1,2,3,4,5,6,7,9,10])))
    negLogl = 1e09;
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal is NaN'])
    keyboard
end

%get prior modes
SetOfpriorModes = SLuniqpair(priorModes);

%(case learnt prior is a von Mises)
%----------------------------------
%if the prior is a von Mises, the unique mode of the prior(s) are those
%found in the databank from (SLMakedatabank).
%In our unimodal experiment it should be [225].
if sum(strcmp(priorShape,'vonMisesPrior'))
    
    %get the modes (should be 225)
    modesPrior1 = SetOfpriorModes(1,:);
    modesPrior2 = SetOfpriorModes(1,:);
    modesPrior3 = SetOfpriorModes(1,:);
    modesPrior4 = SetOfpriorModes(1,:);
    
    %find the trials for each prior
    Prior1 = pstd==80;
    Prior2 = pstd==40;
    Prior3 = pstd==20;
    Prior4 = pstd==10;
    
    %get strength of learnt prior
    klearnt(pstd==80) = fitP(4);
    klearnt(pstd==40) = fitP(5);
    klearnt(pstd==20) = fitP(6);
    klearnt(pstd==10) = fitP(7);
end


%(case learnt prior is bimodal)
%-----------------------------
%get the 4 priors in this order:

if sum(strcmp(priorShape,'bimodalPrior'))
    priorShape = 'bimodalPrior';
    
    %get the modes
    %[145 305]
    %165 285]
    %[185 265]
    %[205 245]
    modesPrior1 = SetOfpriorModes(1,:);
    modesPrior2 = SetOfpriorModes(2,:);
    modesPrior3 = SetOfpriorModes(3,:);
    modesPrior4 = SetOfpriorModes(4,:);
    
    %find the trials for each prior
    Prior1 = priorModes(:,1)==modesPrior1(1,1) & priorModes(:,2)==modesPrior1(1,2);
    Prior2 = priorModes(:,1)==modesPrior2(1,1) & priorModes(:,2)==modesPrior2(1,2);
    Prior3 = priorModes(:,1)==modesPrior3(1,1) & priorModes(:,2)==modesPrior3(1,2);
    Prior4 = priorModes(:,1)==modesPrior4(1,1) & priorModes(:,2)==modesPrior4(1,2);
    
    %The two von Mises that composes the bimodal prior have the same std 20?.
    klearnt(pstd==20) = fitP(6);
end

%-----------
%SET THE LLH
%-----------
%get stimStrengths and strength of evidence
%mean
motdir = unique(d);

%strength
coh24=stimStrength==0.24;
coh12=stimStrength==0.12;
coh06=stimStrength==0.06;
kl(stimStrength==0.24)=fitP(1);
kl(stimStrength==0.12)=fitP(2);
kl(stimStrength==0.06)=fitP(3);
m=1:1:360;

%------------------------------------------------------------
%Bayesian inference with sensory evidence and cardinal priors
%------------------------------------------------------------
%get a matrix (360 possible MAPs, 360 possible motion directions) of
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.
% We only calculate likelihood of the
%data for the actually displayed motion direction and not the full range of
%360 motion directions. It is faster.
%get likelihood of data in condition with stimStrength 24
[~,likeMAP11] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP12] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP13] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP14] = SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 12
[~,likeMAP21] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP22] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP23] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP24] = SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 6
[~,likeMAP31] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior1,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP32] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior2,0,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP33] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior3,0,kcardinal,0,priorShape,TheModel,varargin{:});
[MAP,likeMAP34] = SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior4,0,kcardinal,0,priorShape,TheModel,varargin{:});

%Probability of getting a given estimate given Bayesian inference
%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials)
%for possible values of upo (rows) for each trial (column)
PupoGivenBIwithEvandCard=nan(numel(MAP),numel(d));
for i=1:numel(motdir)
    
    %get displayed motion direction for this trial
    thisd=d==motdir(i);
    
    %get likelihood of data in condition with learnt prior 80
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior1)=likeMAP11(:,motdir(i(ones(sum(thisd&coh24&Prior1),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior1)=likeMAP21(:,motdir(i(ones(sum(thisd&coh12&Prior1),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior1)=likeMAP31(:,motdir(i(ones(sum(thisd&coh06&Prior1),1))));
    
    %get likelihood of data in condition with learnt prior 40
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior2)=likeMAP12(:,motdir(i(ones(sum(thisd&coh24&Prior2),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior2)=likeMAP22(:,motdir(i(ones(sum(thisd&coh12&Prior2),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior2)=likeMAP32(:,motdir(i(ones(sum(thisd&coh06&Prior2),1))));
    
    %get likelihood of data in condition with learnt prior 20
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior3)=likeMAP13(:,motdir(i(ones(sum(thisd&coh24&Prior3),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior3)=likeMAP23(:,motdir(i(ones(sum(thisd&coh12&Prior3),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior3)=likeMAP33(:,motdir(i(ones(sum(thisd&coh06&Prior3),1))));
    
    %get likelihood of data in condition with learnt prior 10
    PupoGivenBIwithEvandCard(:,thisd&coh24&Prior4)=likeMAP14(:,motdir(i(ones(sum(thisd&coh24&Prior4),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh12&Prior4)=likeMAP24(:,motdir(i(ones(sum(thisd&coh12&Prior4),1))));
    PupoGivenBIwithEvandCard(:,thisd&coh06&Prior4)=likeMAP34(:,motdir(i(ones(sum(thisd&coh06&Prior4),1))));
end

%scale to probability.
Z_=sum(PupoGivenBIwithEvandCard);
Z=Z_(ones(size(PupoGivenBIwithEvandCard,1),1),:);
PupoGivenBIwithEvandCard=PupoGivenBIwithEvandCard./Z;

%--------------------------------------------------------
%Bayesian inference with learnt prior and cardinal priors
%--------------------------------------------------------
[~,likeMAP11] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP12] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP13] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP14] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 12
[~,likeMAP21] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP22] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP23] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP24] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

%get likelihood of data in condition with stimStrength 6
[~,likeMAP31] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior1,klearnt1,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP32] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior2,klearnt2,kcardinal,0,priorShape,TheModel,varargin{:});
[~,likeMAP33] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior3,klearnt3,kcardinal,0,priorShape,TheModel,varargin{:});
[MAP,likeMAP34] = SLGirshickBayesLookupTable(ulall,motdir,0,modesPrior4,klearnt4,kcardinal,0,priorShape,TheModel,varargin{:});

%Probability of getting a given estimate
%now get matrix 'PupoGivenBIwithlearntAndCard' of likelihood values
%(upos=1:1:360,trials) for possible values of upo (rows) for each trial
%(column)
PupoGivenBIwithlearntAndCard=nan(numel(MAP),numel(d));
for i=1:numel(motdir)
    
    %get displayed motion direction for this trial
    thisd=d==motdir(i);
    
    %get likelihood of data in condition with learnt prior 80
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior1)=likeMAP11(:,motdir(i(ones(sum(thisd&coh24&Prior1),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior1)=likeMAP21(:,motdir(i(ones(sum(thisd&coh12&Prior1),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior1)=likeMAP31(:,motdir(i(ones(sum(thisd&coh06&Prior1),1))));
    
    %get likelihood of data in condition with learnt prior 40
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior2)=likeMAP12(:,motdir(i(ones(sum(thisd&coh24&Prior2),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior2)=likeMAP22(:,motdir(i(ones(sum(thisd&coh12&Prior2),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior2)=likeMAP32(:,motdir(i(ones(sum(thisd&coh06&Prior2),1))));
    
    %get likelihood of data in condition with learnt prior 20
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior3)=likeMAP13(:,motdir(i(ones(sum(thisd&coh24&Prior3),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior3)=likeMAP23(:,motdir(i(ones(sum(thisd&coh12&Prior3),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior3)=likeMAP33(:,motdir(i(ones(sum(thisd&coh06&Prior3),1))));
    
    %get likelihood of data in condition with learnt prior 10
    PupoGivenBIwithlearntAndCard(:,thisd&coh24&Prior4)=likeMAP14(:,motdir(i(ones(sum(thisd&coh24&Prior4),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh12&Prior4)=likeMAP24(:,motdir(i(ones(sum(thisd&coh12&Prior4),1))));
    PupoGivenBIwithlearntAndCard(:,thisd&coh06&Prior4)=likeMAP34(:,motdir(i(ones(sum(thisd&coh06&Prior4),1))));
end

%scale to probability.
Z_=sum(PupoGivenBIwithlearntAndCard);
Z=Z_(ones(size(PupoGivenBIwithlearntAndCard,1),1),:);
PupoGivenBIwithlearntAndCard=PupoGivenBIwithlearntAndCard./Z;

%get probabilities of the range of percepts "upo" given random estimation
%------------------------------------------------------------------------
PupoGivenRand=ones(360,numel(d))/360;

%calculate the fraction of trials that is controlled by each process
%We use a divisive normalization competition rule between the likelihood,
%the learnt prior and the cardinal prior. Basically the stronger the
%representation and the
%most often it will be chosen. We use the strength of the von Mises
%composing the cardinal prior and not some measure of the overall strength
%of the cardinal prior. But both are correlated. For fixed modes, if the
%strength of the von Mises (one value here) increases, so does the strength
%of the overall distribution. I think...
%The probability of random estimation is fixed.
weightPriorCardmn   = klearnt./(klearnt+kl);
weightLlhCardmn     = kl./(klearnt+kl);

%scale the mixing weights to probabilities (all sum to 1)
%this is not equal to 1. We want to make it equal to 1.
sumP=weightPriorCardmn + weightLlhCardmn + Prandom;
PpriorCardmn     = 1 - (Prandom + weightLlhCardmn)./sumP;
PllhCardmn       = 1 - (Prandom + weightPriorCardmn)./sumP;
Prandnew     = unique(1 - (weightLlhCardmn + weightPriorCardmn)./sumP);

%take just one value of Prandnew because sometimes because of numerical
%instability there is slight variation in the value of Prandnew which
%produces more than one value. Mathematically there should be only one
%value because priors and llh strength have been scaled some that they sum
%to 1 and Prandom is a constant. Thus Prandnew must be one value.
Prandnew=Prandnew(1);

%repeat P(choose prior), matrix values, of each trial (columns) for each
%possible data value (rows).
numm = numel(m);
PpriorCardmnall=PpriorCardmn(ones(numm,1),:);
PllhCardmnall=PllhCardmn(ones(numm,1),:);

%We assume here a particular model of competition (llh mean vs prior mean).
%Either one chooses evidence (llh mean) or he chooses the prior mean.
%P(choose evidence|motion) = P(obs evidence & choose evidence/motion and prior)
%P(choose prior mean|motion) = P(choose prior mean/motion and prior)
%Thus equal llh and prior strengths do not produces a bimodal estimate
%density with two peaks with same amplitude. It would be the case if either
%subject chose the motion direction (mean of sensory density) or the prior
%mean(Competition between prior and sensory density). Here it's competition
%between likelihood mean and prior mean. An alternative model with make
%motion direction (mean of evidence density) and prior mean compete instead
%of the mean of likelihood (sampled at each trial from evidence density)
%and the prior mean.
PupoGivenModel = PupoGivenBIwithEvandCard.*PllhCardmnall + PupoGivenBIwithlearntAndCard.*PpriorCardmnall + PupoGivenRand.*Prandnew;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel))==1
    keyboard
    fprintf('%s \n','Something s wrong. PupoGivenModel are probabilties and should sum to 1')
end

%convolve with motor noise
%--------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo=0:1:359;
Pmot=vmPdfs(upo,0,km,'norm');
Pmotforconv=Pmot(:,ones(1,numel(d)));
PestimateGivenModel=SLcircConv(PupoGivenModel,Pmotforconv);

%get read of extremely small negative numbers that occurs because of
%instability numerical instability. Set them to zero.
%PestimateGivenModel(PestimateGivenModel<0)=0;

%check that probability of estimates Given Model are positive values.
%circular convolution sometimes produces negative values very close to zero
%(order of -10^-18). Those values produce infinite -log likelihood which are
%only need a single error in estimation to be rejected during model fitting
%(one trial has +inf -log likelihood to be predicted by the model).
%This is obviously too conservative. We need a minimal non-zero lapse rate.
%Prandom that allows for errors in estimation. So we add 10^320 the minimal
%natural number available in matlab. This means that every time an estimate
%that is never produced by the model wfithout lapse rate is encountered
%-loglikelihood increases by -log(10^-320) = 737 and is thus less likely to
%be a good model (lowest -logLLH). But the model is not rejected altogether
%(it would be -log(0) = inf). In the end models that cannot account for
%error in estimates are more likely to be rejected than models who can.
if ~isempty(find(PestimateGivenModel<=0))
    
    PestimateGivenModel(PestimateGivenModel<=0) = 10^-320;
    fprintf('(SLgetLoglBayesianModel) P(estimates|model) were sometimes negative but very close to zero < 10^-10 and were thus floored at the lowest positive numerical value 10^-320. \n')
    
end

%set upo to initial values any case we use it later
upo=1:1:360;

%trial-predictions (sample estimate density). Variability in estimate
%should reflect variability in measurement density and motor noise
for i=1:length(d)
    pred(i)=randsample(1:1:360,1,'true',PestimateGivenModel(:,i));
end

%draw sampled predictions
%------------------------
%draw predictions and data
%drawCircStat(pred,d,coh,pstd);
%drawCircStat(data,d,coh,pstd);

%get loglikelihood of data
%-------------------------
%single trial's measurement, its position(row) for each trial(col) and its
%probability (also maxlikelihood of trial's data). Checked many times. It
%works.
%make sure sub2ind inputs are the same size
if sum(size(data)~=size([1:1:numel(d)]'))==2
    data=data';
end
idx = sub2ind(size(PestimateGivenModel),data,[1:1:numel(d)]');
PdataGivenModel = PestimateGivenModel(idx);

%We use log likelihood because likelihood is so small that matlab cannot
%encode it properly (numerical unstability). We can use single trials log
%likelihood to calculate AIC in the conditions that maximize differences in
%predictions of two models.
Logl_pertrial = log(PdataGivenModel);

%We use -sum(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
%logL = -sum(logL_pertrial);

%We use -sum(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
%negSumlogL=-sum(log(PdataGivenModel));

%We use E(log likelihood) as an objective function to minimize with
%matlab non linear optimization search.
negLogl = -sum(Logl_pertrial);

%Look at fitting. It is 3X faster without drawing.
ti = toc(ticLogL);
if sum(strcmp(varargin{1},'vonMisesPrior')) && sum(strcmp(varargin{1},'bimodalPrior'))
else
    fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
        negLogl,kl1,kl2,kl3,klearnt1,klearnt2,klearnt3,klearnt4,kcardinal,Prandom,km,ti)
end

%logl for combined priors
function [negLogl,fitP,Logl_pertrial] = getLoglCombinedPriors(dataVM,dataBim,...
    displVM,displBim,stimStrengthVM,stimStrengthBim,pstdVM,pstdBim,fitP,...
    priorModesVM,priorModesBim,TheModel,...
    varargin)

%we use the same fit parameters to get log likelihood for von Mises prior
%data and bimodal prior data separately. We then calculated the logl of the
%overall data set by summing the two.

ticLogL = tic;

%negLogl
%von Mises prior
[negLoglVM,fitP,Logl_pertrialVM] = getLogl(dataVM,displVM,stimStrengthVM,pstdVM,...
    fitP,'vonMisesPrior',priorModesVM,TheModel,varargin{1});

if isempty(Logl_pertrialVM)
    displ('empty')
end

%bimodal prior
fitP(4:7) = fitP(6);
[negLoglBim,fitP,Logl_pertrialBim] = getLogl(dataBim,displBim,stimStrengthBim,...
    pstdBim,fitP,'bimodalPrior',priorModesBim,TheModel,varargin{1});

%combine
negLogl = sum([negLoglVM;negLoglBim]);
Logl_pertrial = [Logl_pertrialVM; Logl_pertrialBim];


%Look at fitting.
ti = toc(ticLogL);
fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
    negLogl,fitP(1),fitP(2),fitP(3),fitP(4),fitP(5),fitP(6),fitP(7),fitP(8),fitP(9),fitP(10),ti)



%-----------------------------------------------------------
%Case cross-validated R-squared fit (estimates mean and std)
%-----------------------------------------------------------
%valid only for lots of data
function [pred,data,disp,stimStrength,pstd,stdPa,output] = CrossValR2fit(databank,...
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
fprintf('%s \n',['(CrossValR2fit) Now fitting model to data mean and std',...
    ' for each experimental condition'])

%(case check fitting fake data). Get simulated data in the workspace.
if sum(strcmp(varargin{1},'testFakedata'))
    data=evalin('base','pred');
    disp=evalin('base','d');
    stimStrength=evalin('base','stimStrength');
    pstd=evalin('base','pstd');
    priorModes=evalin('base','priorModes');
    
else
    
    %subjects' data
    data=round(cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1)));
    
    %conditions
    disp=cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
    stimStrength=cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd=cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes=cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    data(data==0)=360;
    
    %remove missing data
    data=data(isnan(data)==0);
    disp=disp(isnan(data)==0);
    stimStrength=stimStrength(isnan(data)==0);
    pstd=pstd(isnan(data)==0);
end

%store data and conditions
output.data = data;
output.disp = disp;
output.pstd = pstd;
output.priorModes = priorModes;

%----------------------------------
%1 - Sort data as training and test
%----------------------------------
%Divide data into 2 sets, gets best fit parameters on 1 sets (training) and
%use them to calculate the R2 of the model fitted to the test data.
%We chose 2 sets and not more because the test set must contain enough
%samples at the tails to be able to compute std. e,g., when there is
%only one sample trial std =0 and is not a proper estimate of the estimate.
%estimate fall down on the plot at the tail. This artefact impairs more
%competition model at the tail than Sampling model because Competition
%model raises more std at the tails than Sampling model which is what
%the data are doing when we plot the entire dataset.
fprintf('%s \n','(CrossValR2fit) Sorting training and test sets ...')
numSet = 5;
output.CrossVR2 = nan(numSet,1);

%loop over sets
for thisSet = 1 : numSet
    
    %status
    fprintf('%s \n',['(CrossValR2fit) Calculating R-squared for set ',...
        num2str(thisSet),'/',num2str(numSet)])
    
    %------------------------------------------------------------
    %1 - Divide data into 2 sets that each contain all conditions
    %------------------------------------------------------------
    %(case von Mises)
    %----------------
    if strcmp(priorShape,'vonMisesPrior')
        
        %each trial's experimental condition
        [~,~,idxCondAlltrials] = SLuniqpair([pstd stimStrength disp]);
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
        stimStrengthTrain = stimStrength(TrialsTrain);
        pstdTrain = pstd(TrialsTrain);
        priorModesTrain = priorModes(TrialsTrain);
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
        
        %training data mean and std are sorted by conditions
        [meanData,stdData,myF] = SLCircStat(dataTrain,dispTrain,...
            stimStrengthTrain,pstdTrain);
        meanData = meanData(:);
        stdData = stdData(:);
        myF1 = myF.f1.D3(:);
        myF2 = myF.f2.D3(:);
        myF3 = myF.f3.D3(:);
        myCond = [myF3 myF2 myF1];
        
        %data descriptive stats
        statsTrain = SLmakeCircStat(dataTrain,pstdTrain,stimStrengthTrain,dispTrain);
        
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
        fprintf('%s \n',['(CrossValR2fit) Initialize your ',...
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
        options = optimset('MaxIter',output.numIter,'MaxFunEvals',output.numFunEval);
        fprintf('%s \n',['(CrossValR2fit) Fitting with ',num2str(output.numIter),...
            ' iterations and ',num2str(output.numFunEval),' function evaluations...'])
        
        %loop over 10 sets of initial parameters
        vrg = varargin{1};
        %parfor i = 1 : size(k0,1)
        for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n',['(CrossValR2fit)',' Set of initial parmeters: ',...
                num2str(i),'/',num2str(size(k0,1))])
            
            %Nelder-Mead (simplex search algorithm)
            [fitPtmp,SSE,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                makeSSE(dispTrain,...
                stimStrengthTrain,...
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
        fprintf('%s \n',['(CrossValR2fit) Identifying minimum SSE and',...
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
            output.fitP.nm={'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10',...
                'kcardinal','Prand','km'};
        end
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'bimodalPrior')
            output.fitP.nm={'coh24','coh12','coh6','pstd145_305','pstd165_285',...
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
        
        SSE = makeSSE(dispTrain,stimStrengthTrain,pstdTrain,fitPinput,...
            priorShape,priorModesTrain,TheModel,meanData,stdData,...
            dataCond,varargin{1});
        
        %SSE
        output.minSSE = SSE;
        
        %fitP
        output.fitPbkp{thisSet} = fitPinput;
        output.fitP.p(thisSet,:) = fitPinput;
        output.fitP.nm = {'coh24','coh12','coh6','klearnt80','klearnt40','klearnt20',...
            'klearnt10','Kcardinal','Prand','km'};
        pred=[];
        stdPa=[];
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
    cohTest  = stimStrength(output.TrialsTest{thisSet});
    pstdTest = pstd(output.TrialsTest{thisSet});
    priorModesTest = priorModes(output.TrialsTest{thisSet});
    
    %mean and std of test data
    [meanDataTest,stdDataTest,myFTest] = SLCircStat(dataTest,dispTest,...
        cohTest,pstdTest);
    meanDataTest = meanDataTest(:);
    stdDataTest = stdDataTest(:);
    myF1 = myFTest.f1.D3(:);
    myF2 = myFTest.f2.D3(:);
    myF3 = myFTest.f3.D3(:);
    myCondTest = [myF3 myF2 myF1];
    
    %data descriptive stats
    statsTest = SLmakeCircStat(dataTest,pstdTest,cohTest,dispTest);
    
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
    drawMeanPreCentered(meanDataTest,stdDataTest,[],[],myCondTest)
    
    %test SSE
    SSE = makeSSE(dispTest,cohTest,pstdTest,output.fitP.p(thisSet,:),priorShape,...
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

output.CrossVmeanR2 = mean(output.CrossVR2);

%-----------------------------------------------
%Case cross-validated R-squared fit (trial-data)
%-----------------------------------------------
function [pred,data,disp,stimStrength,pstd,stdPa,output] = CrossValR2FitTrialData(...
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
    stimStrength = evalin('base','stimStrength');
    pstd = evalin('base','pstd');
    priorModes = evalin('base','priorModes');
    
else
    
    %subjects' data
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1)));
    
    %conditions
    disp = cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
    stimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd  =cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    stimStrength = stimStrength(isnan(data)==0);
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
numSet = 5;
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
        [~,~,idxCondAlltrials] = SLuniqpair([pstd stimStrength disp]);
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
        cohTrain = stimStrength(TrialsTrain);
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
            cohTrain,pstdTrain);
        meanData = meanData(:);
        stdData = stdData(:);
        myF1 = myF.f1.D3(:);
        myF2 = myF.f2.D3(:);
        myF3 = myF.f3.D3(:);
        myCond = [myF3 myF2 myF1];
        
        %data descriptive stats
        %statsTrain = SLmakeCircStat(dataTrain,pstdTrain,cohTrain,dispTrain);
        
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
        fprintf('%s \n',['(CrossValR2fit) Fitting with ',...
            num2str(output.numIter),...
            ' iterations and ',num2str(output.numFunEval),...
            ' function evaluations...'])
        
        %loop over 10 sets of initial parameters
        vrg = varargin{1};
        
        %parfor i = 1 : size(k0,1)
        for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n',['(CrossValR2fit)',' Set of initial parmeters: ',...
                num2str(i),'/',num2str(size(k0,1))])
            
            %Nelder-Mead (simplex search algorithm)
            [fitPtmp,SSE,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                makeSSE(dispTrain,...
                cohTrain,...
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
            output.fitP.nm={'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10',...
                'kcardinal','Prand','km'};
        end
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'bimodalPrior')
            output.fitP.nm={'coh24','coh12','coh6','pstd145_305','pstd165_285',...
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
        SSE = makeSSE(dispTrain,cohTrain,pstdTrain,fitPinput,...
            priorShape,priorModesTrain,TheModel,'Trial',dataTrain,meanData,...
            stdData,dataCond,output,varargin{1});
        
        %SSE
        output.minSSE = SSE;
        
        %fitP
        output.fitPbkp{thisSet} = fitPinput;
        output.fitP.p(thisSet,:) = fitPinput;
        output.fitP.nm = {'coh24','coh12','coh6','klearnt80','klearnt40','klearnt20',...
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
    cohTest  = stimStrength(output.TrialsTest{thisSet});
    pstdTest = pstd(output.TrialsTest{thisSet});
    priorModesTest = priorModes(output.TrialsTest{thisSet});
    
    %mean and std of test data
    [meanDataTest,stdDataTest,myFTest] = SLCircStat(dataTest,dispTest,...
        cohTest,pstdTest);
    meanDataTest = meanDataTest(:);
    stdDataTest = stdDataTest(:);
    myF1 = myFTest.f1.D3(:);
    myF2 = myFTest.f2.D3(:);
    myF3 = myFTest.f3.D3(:);
    myCondTest = [myF3 myF2 myF1];
    drawMeanPreCentered(meanDataTest,stdDataTest,[],[],myCondTest)
    
    
    %data descriptive stats
    statsTest = SLmakeCircStat(dataTest,pstdTest,cohTest,dispTest);
    
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
    drawMeanPreCentered(meanDataTest,stdDataTest,[],[],myCondTest)
    
    %test SSE
    SSE = makeSSE(dispTest,cohTest,pstdTest,output.fitP.p(thisSet,:),priorShape,...
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

%-------------------------------------------------------
%Case cross-validated max loglikelihood fit (trial-data)
%-------------------------------------------------------
function [pred,data,disp,stimStrength,pstd,stdPa,output] = CrossValR2FitTrialDataMaxLL(...
    databank,...
    fitPinput,...
    initp,...
    priorShape,...
    filename,...
    output,...
    varargin)

%STEPS of the procedure
%1. First, sort data as training and test.

%2. Second, get the best fit parameters that minimize - loglikelihood (i.e,?
% the maximum loglikelihood between training data and predictions.

%3. Use the fit Parameters to generate predictions for the test data set
%and calculate test loglikelihood.

%4. Start a new training data set and repeat the procedure.

%5. Store the R2 for each of the 5 sets and calculate the mean loglikelihood.

%status
fprintf('%s \n','(CrossValR2FitTrialDataMaxLL) Now fitting model to trial data')

%(case check fitting fake data).
%Get simulated data in the workspace.
if sum(strcmp(varargin{1},'testFakedata'))
    
    data = evalin('base','pred');
    disp = evalin('base','d');
    stimStrength = evalin('base','stimStrength');
    pstd = evalin('base','pstd');
    priorModes = evalin('base','priorModes');
    
else
    
    %subjects' data
    data = round(cell2mat(databank.data(:,(strcmp(databank.nm,'est_dir'))==1)));
    
    %conditions
    disp = cell2mat(databank.data(:,(strcmp(databank.nm,'sample_dir'))==1));
    stimStrength = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
    pstd  =cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
    priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));
    
    %make sure 360 and 0 degrees are same
    data(data==0) = 360;
    
    %remove missing data
    data = data(isnan(data)==0);
    disp = disp(isnan(data)==0);
    stimStrength = stimStrength(isnan(data)==0);
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
%2 sets is best when we don't have enough data to have non-noisy test set.
fprintf('%s \n','(CrossValR2FitTrialDataMaxLL) Sorting training and test sets ...')
numSet = 2;
output.CrossVR2 = nan(numSet,1);

%loop over sets
for thisSet = 1 : numSet
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Calculating -logl for set ',...
        num2str(thisSet),'/',num2str(numSet)])
    
    %------------------------------------------------------------
    %1 - Divide data into 2 sets that each contain all conditions
    %------------------------------------------------------------
    %(case von Mises)
    %----------------
    if strcmp(priorShape,'vonMisesPrior')
        
        %each trial's experimental condition
        [~,~,idxCondAlltrials] = SLuniqpair([pstd stimStrength disp]);
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
        cohTrain = stimStrength(TrialsTrain);
        pstdTrain = pstd(TrialsTrain);
        priorModesTrain = priorModes(TrialsTrain);
    end
    
    %store
    output.dataTrainBkp{thisSet} = dataTrain;
    
    %------------------------------------------------
    %2 - max likelihood fit of model to training data
    %-------------------------------------------------
    %(case von Mises)
    %----------------
    %Calculate data mean and std
    if strcmp(priorShape,'vonMisesPrior')
        
        %status
        fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Calculating training data means'...
            ' and stds per condition for set ',...
            num2str(thisSet),'/',num2str(numSet)])
        
        %training data mean and std are sorted by conditions
        [meanData,stdData,myF] = SLCircStat(dataTrain,dispTrain,...
            cohTrain,pstdTrain);
        meanData = meanData(:);
        stdData = stdData(:);
        myF1 = myF.f1.D3(:);
        myF2 = myF.f2.D3(:);
        myF3 = myF.f3.D3(:);
        myCond = [myF3 myF2 myF1];
        
        %data descriptive stats
        %statsTrain = SLmakeCircStat(dataTrain,pstdTrain,cohTrain,dispTrain);
        
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
        k0(1,:) = [kllh    klearnt      kcardinal fractRand motorNoise];
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
        k0(10,:) = [repmat(nanmean(k0(1,1:8)),1,7) kcardinal fractRand motorNoise];
        
        %status
        output.numInitP = size(k0,1);
        fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Initialize your ',...
            num2str(output.numInitP),' initial parameters'])
        
        %(case we do not want to fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            TheModel = 'withoutCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) The model does not contain',...
                'a cardinal prior'])
        else
            TheModel = 'withCardinal';
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) The model contains',...
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
        
        %fitting options
        %--------------
        %for default fitting
        output.numIter = [];
        output.numFunEval = [];
        %options = [];
        
        %for fast debugging or fitting
        %output.numIter = 100;
        %output.numFunEval = 100;
        %options = optimset('MaxIter',output.numIter,'MaxFunEvals',...
        %    output.numFunEval);
        
        %for debugging
        options = optimset('MaxIter',1,'MaxFunEvals',1);
        
        %status
        fprintf('%s \n',['(CrossValR2fit) Fitting with ',...
            num2str(output.numIter),...
            ' iterations and ',num2str(output.numFunEval),...
            ' function evaluations...'])
        
        %loop over 10 sets of initial parameters
        vrg = varargin{1};
        
        %parfor i = 1 : size(k0,1)
        for i = 1 : size(k0,1)
            
            %status
            fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL)',' Set of initial parmeters: ',...
                num2str(i),'/',num2str(size(k0,1))])
            
            %Nelder-Mead (simplex search algorithm)
            [fitPtmp,negLogl,exitflag,outputFit] = fminsearch(@(fitPtmp) ...
                getLogl(dataTrain,...
                dispTrain,...
                cohTrain,...
                pstdTrain,...
                fitPtmp,...
                priorShape,...
                priorModesTrain,...
                TheModel,vrg),...
                k0(i,:),...
                options);
            
            %Fit parameters and SSE
            fitPbkp(i,:)    = fitPtmp;
            negLogl_bkp(i)  = negLogl;
            outputFitBkp(i) = {outputFit};
            exitflagBkp{i}  = exitflag;
        end
        
        %backups
        output.fitPbkp{thisSet}     = fitPbkp;
        output.negLogl_bkp{thisSet} = negLogl_bkp;
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
        
        [minNegLogL,position] = min(negLogl_bkp(:));
        i = ind2sub(size(negLogl_bkp),position);
        
        %minimal SSE
        output.minNegLogL = minNegLogL;
        
        %optimal parameters
        output.fitP.p(thisSet,:) = fitPbkp(i,:);
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'vonMisesPrior')
            output.fitP.nm={'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10',...
                'kcardinal','Prand','km'};
        end
        
        %(case von Mise prior)
        %---------------------
        if strcmp(priorShape,'bimodalPrior')
            output.fitP.nm={'coh24','coh12','coh6','pstd145_305','pstd165_285',...
                'pstd185_265','pstd205_245',...
                'kcardinal','Prand','km'};
        end
        
        %(case we did not fit a cardinal prior). Kcardinal=NaN.
        if isnan(kcardinal)
            output.fitP.p(thisSet,strcmp(output.fitP.nm,'kcardinal'))=NaN;
        end
        
        %predictions
        pred=[];
        
        %std of model parameters (interior-point).
        stdPa = [];
        
        %save data in case bug
        save(filename)
    end
    
    %----------------------------------------------------------------------
    %2 -  Use training best fit parameters to calculate R2 in test data set
    %----------------------------------------------------------------------
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Use training best fit parameters to',...
        'Calculate max loglikelihood in test data...'])
    
    %make test data set
    %select test trials (5th set left out for test)
    output.TrialsTest{thisSet} = output.CrossVSet==testSetID;
    dataTest = data(output.TrialsTest{thisSet});
    output.dataTestBkp = dataTest;
    
    %test conditions
    dispTest = disp(output.TrialsTest{thisSet});
    cohTest  = coh(output.TrialsTest{thisSet});
    pstdTest = pstd(output.TrialsTest{thisSet});
    priorModesTest = priorModes(output.TrialsTest{thisSet});
    
    %mean and std of test data
    [meanDataTest,stdDataTest,myFTest] = SLCircStat(dataTest,dispTest,...
        cohTest,pstdTest);
    meanDataTest = meanDataTest(:);
    stdDataTest = stdDataTest(:);
    myF1 = myFTest.f1.D3(:);
    myF2 = myFTest.f2.D3(:);
    myF3 = myFTest.f3.D3(:);
    myCondTest = [myF3 myF2 myF1];
    
    %draw test data
    drawMeanPreCentered(meanDataTest,stdDataTest,[],[],myCondTest)
    
    %data descriptive stats
    statsTest = SLmakeCircStat(dataTest,pstdTest,cohTest,dispTest);
    
    %status
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) IMPORTANT !!!! Please Check',...
        ' the test data that is plotted. e.g.,number of sample data for',...
        ' a condition may not be enough. Not enough data will provide',...
        ' poor R-squared !'])
    
    %display descriptive stats to check before fitting
    %M = [statsTest.count statsTest.conditions];
    %fprintf('%8.0f %8.0f %8.2f %8.0f \n',M')
    
    %back up
    output.statsTest{thisSet} = statsTest;
    
    %ckeck data mean and std visually
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) IMPORTANT !!!! Now drawing data',...
        ' from test set....'])
    drawMeanPreCentered(meanDataTest,stdDataTest,[],[],myCondTest)
    
    %test max loglikelihood
    output.CrossVnegLogl(thisSet) = getLogl(dataTest,...
        dispTest,...
        cohTest,...
        pstdTest,...
        output.fitP.p(thisSet,:),...
        priorShape,...
        priorModesTest,...
        TheModel,varargin{1});
    output.CrossVmaxLogl(thisSet) = - output.CrossVnegLogl(thisSet);
    
    %status
    fprintf('%s \n','----------------------------------------------------')
    fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) maxLogl for set ',...
        num2str(thisSet),'/',num2str(numSet),'= ',...
        num2str(output.CrossVmaxLogl(thisSet))])
    fprintf('%s \n','----------------------------------------------------')
end

%mean max loglikelihood
%status
fprintf('%s \n',['(CrossValR2FitTrialDataMaxLL) Calculating mean logl over the',...
    num2str(numSet),' test sets'])

output.CrossVmeanMaxLogl = mean(output.CrossVmaxLogl(thisSet));

%SSE (trial and average)
function [SSE,fitP,output] = makeSSE(displ,coh,pstd,fitP,priorShape,priorModes,...
    TheModel,TrialOrMean,trialData,meanData,stdData,dataCond,output,varargin)

%time
ticSSE = tic;

%remove missing conditions (NaN)
pos =~ isnan(meanData);
meanData = meanData(pos);
stdData  = stdData(pos);
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
if fitP(9)>1
    SSE = 1e9;
    return
end
if any(fitP<0)
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
[meanPred,stdPred,cond,~,~,output] = SLmakePredictionsCompetitionModel(displ,stimStrength,pstd,fitP,...
    priorShape,priorModes,TrialOrMean,output,varargin{:});


%case least-square fit to estimate mean and std
%----------------------------------------------
if strcmp(TrialOrMean,'Mean')
    
    %check prediction and data condition do not match
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
    
    %(Visual check) estimate mean and std and predictions
    drawMeanPreCentered(meanData,stdData,meanPred,stdPred,dataCond)
end


%case least-square fit to trial data
%-----------------------------------
if strcmp(TrialOrMean,'Trial')
    trialPred = output.TrialPred;
    Error = SLvectors2signedAngle(trialData,trialPred);
    
    %%(Visual check) data mean and pred
    %%data
    %[meanData,stdData,myF] = SLCircStat(trialData,displ,...
    %    coh,pstd);
    %meanData = meanData(:);
    %stdData = stdData(:);
    %myCond = [myF.f3.D3(:) myF.f2.D3(:) myF.f1.D3(:)];
    
    %predictions
    %predictions of std are not possible because model outputs the same mean
    %for each same condition.
    %meanPred = SLCircStat(trialPred,displ,...
    %    coh,pstd);
    %meanPred = meanPred(:);
    
    %draw
    %drawMeanPreCentered(meanData,stdData,meanPred,[],myCond)
end

%SSE
SSE = nansum(Error.^2);

%debug if SSE is a complex value
if ~isreal(SSE)
    dbstop
end

%fit info
ti = toc(ticSSE);
fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f  %.2f \n',...
    SSE,fitP,ti)

%R^2 (trial and average)
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
%data mean and std for each experimental condition
function [meanData,stdData,dataCond] = makeDataMeanAndStd(data,disp,...
    stimStrength,pstd,priormodes,priorShape)

%(case von Mises prior)
%----------------------
%Calculate data mean and std for each experimental condition (disp,stimStrength,pstd)
if strcmp(priorShape,'vonMisesPrior')
    
    %data are sorted by experimental conditions
    [meanData,stdData,myF] = SLCircStat(data,disp,stimStrength,pstd);
    meanData = meanData(:);
    stdData = stdData(:);
    myF1 = myF.f1.D3(:);
    myF2 = myF.f2.D3(:);
    myF3 = myF.f3.D3(:);
    dataCond = [myF3 myF2 myF1];
end

%remove missing conditions (NaN) in the data
pos =~ isnan(meanData);
meanData = meanData(pos);
stdData = stdData(pos);
dataCond = dataCond(pos,:);

%data distributions per condition
function [p,xpdf] = makeDataDist(data,d,stimStrength,pstd,priorModes,cond,priorShape)

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
for i = 1 : size(cond,1)
    
    %this condition
    %case von Mises prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        thisCon = pstd==cond(i,1) & stimStrength==cond(i,2) & d==cond(i,3);
        
        %case bimodal prior
        %--------------------
    elseif strcmp(priorShape,'bimodalPrior')
        thisCon = priorCond==cond(i,1) & stimStrength==cond(i,2) & d==cond(i,3);
        
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
    
    %show sample size
    fprintf('%s \n', ['(analyses) ',num2str(numel(dtoHist)),' trials'])
end




%----
%draw
%----
%draw model predictions mean and std
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
F.f2.nm = 'stimStrength';
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

%draw model predictions distribution
function drawDataAndPreDist(Pdata,bins,Ppred,d,stimStrength,pstd,priorModes,...
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
Thecohs = unique(stimStrength);
dlinDisttoPrior = d - priorModes;

%number of conditions
numPriors = numel(Thepriors);
numCoh = numel(Thecohs);

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

%initialize the conditions to plot in each figure (coh)
condThisCoh = nan(numDir*numPriors,size(cond,2));
condThisCoh(:,1) = SLreplicateRows(Thepriors,numDir);
condThisCoh(:,3) = repmat(motDir,4,1);
numAllCond = size(condThisCoh,1);

%this stimStrength
for j = 1 : numCoh
    
    %figure
    figure('color','w',...
        'Position',[1 scrsz(4)/2 scrsz(3)/width scrsz(4)])
    
    %this stimStrength
    thiscoh = Thecohs(j);
    
    %sync subplots and cond
    condThisCoh(:,2) = repmat(Thecohs(j),numDir*numPriors,1);
    for i = 1 : numAllCond
        
        %time
        t1 = tic;
        
        %axis
        hs(i) = subplot(numDir,numPriors,dirpos(i));
        
        %condition
        thisCon = SLfindRow(condThisCoh(i,:),cond);
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
            
            title(['(',num2str(thiscoh*100),'% coh -',...
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

%draw predicted mean and std (centered at prior mean)
function drawMeanPreCentered(meanData,stdData,meanPred,stdPred,dataCond)

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
F.f2.nm = 'stimStrength';
F.f2.L = unique(F.f2.i);
F.f2.L = sort(F.f2.L,'descend');
F.f2.n = numel(F.f2.L);

F.f3.i = dataCond(:,1);
F.f3.nm = 'Prior std';
F.f3.L = unique(F.f3.i);
F.f3.L = sort(F.f3.L,'descend');
F.f3.n = numel(F.f3.L);

%Graphics
F.f2.color = {[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};

F.f2.colorPre = {[0.2 0 0],...
    [0.97 0.2 0],...
    [0.8 0.4 0],...
    [0.3 0.3 0]};

%Mean data and predictions
scrsz = get(0,'ScreenSize');
fig1 = figure(1);
set(fig1,'color','w','Position',[0 1.4*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
ymax = nan(F.f2.n,1);
ymin = nan(F.f2.n,1);
h = nan(F.f2.n,1);
count = 0;

%graphics
marksz = 100;
ftsz = 12;

for j = 1 : F.f2.n
    
    h(j) = subplot(1,F.f2.n,j);
    axis square
    
    for i = 1 : F.f3.n
        
        %data for this conditions
        thisC = F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %recentering to prior
        %--------------------
        %recenter all data relative to prior
        %calculate distance to prior mean 225 deg
        %sort distance
        xCentered = SLvectors2signedAngle(F.f1.i( thisC ),225);
        [xCenteredSorted,IA] = sort(xCentered,'ascend');
        
        %predictions
        yCentered = SLvectors2signedAngle(meanPred( thisC ),225);
        %yCentered = meanPred( thisC );
        yCenteredSorted = yCentered(IA);
        
        %data
        yDataCentered = SLvectors2signedAngle(meanData( thisC ),225);
        %yDataCentered = meanData( thisC );
        yDataCenteredSorted = yDataCentered(IA);
        
        %sort positions and labels for plot such that 225 deg the prior
        %mean is at the center of the plot (when distance = 0).
        %Sorting is done both for x and y data such that they remain
        %aligned across conditions.
        xTickCentered = SLvectors2signedAngle(F.f1.L,225);
        [xTickCenteredSorted,I] = sort(xTickCentered,'ascend');
        yTickCenteredSorted = xTickCenteredSorted;
        xtickLabel = F.f1.L(I);
        ytickLabel = xtickLabel;
        
        %plot
        count = count+1;
        myPlot1(count) = scatter(xCenteredSorted,yDataCenteredSorted,...
            marksz,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.f2.color{i},...
            'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));
        
        myPlot2(count) = plot(xCenteredSorted,yCenteredSorted,...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
        
        %graphics
        ymax(j,i) = max([yDataCenteredSorted; yCenteredSorted]);
        ymin(j,i) = min([yDataCenteredSorted; yCenteredSorted]);
    end
    
    %x and ylabel
    if j==1
        ylabel('Predicted average estimate (deg)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
    
    %centered x labels
    set(gca,'ytick',yTickCenteredSorted(3:8:end),...
        'yticklabel',ytickLabel(3:8:end))
    %set(gca,'ytick',F.f1.L(7:8:end),...
    %    'yticklabel',F.f1.L(7:8:end))
    set(gca,'xtick',xTickCenteredSorted(3:8:end),...
        'xticklabel',xtickLabel(3:8:end))
    
    %selected legend
    lg = legend(myPlot2(j));
    legend(lg,'boxoff')
end

%x and ylimits
set(h,'ylim',[min(ymin(:)) max(ymax(:))])
set(h,'xlim',[min(ymin(:)) max(ymax(:))])

%clear up
SLremoveDeadSpace;
%SLConventionUp(gcf)

%----
%std
%----
fig2 = figure(2);
set(fig2,'color','w','Position',[0 0 scrsz(3)/3 scrsz(4)/3])
ymax = nan(F.f2.n,1);
ymin = nan(F.f2.n,1);
h = nan(F.f2.n,1);
count = 0;
for j = 1 : F.f2.n
    
    h(j) = subplot(1,F.f2.n,j);
    axis square
    
    for i = 1 : F.f3.n
        
        %this conditions data
        thisC = F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %transform to distance relative to prior
        %this automatically align prior mean to the center of the plot
        %(distance=0)
        %distance to prior mean
        %sort distance to prior mean
        xCentered = SLvectors2signedAngle(F.f1.i( thisC ),225);
        [xCenteredSorted,IA] = sort(xCentered,'ascend');
        
        %data
        yDataCentered = stdData( thisC );
        yDataCenteredSorted = yDataCentered(IA);
        
        %predictions
        yCentered = stdPred( thisC );
        yCenteredSorted = yCentered(IA);
        
        %sort positions and labels for plot
        xTickCentered = SLvectors2signedAngle(F.f1.L,225);
        [xTickCenteredSorted,I] = sort(xTickCentered,'ascend');
        xtickLabel = F.f1.L(I);
        
        %plot
        count = count+1;
        
        myPlot1(count) = scatter(xCenteredSorted,yDataCenteredSorted,...
            marksz,...
            'MarkerEdgeColor','w',...
            'MarkerFaceColor',F.f2.color{i},...
            'displayname',strcat(F.f3.nm,':',num2str(F.f3.L(i))));
        
        
        myPlot2(count) = plot(xCenteredSorted,yCenteredSorted,...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
        
        %graphics
        xmax(j,i) = max(xCenteredSorted);
        xmin(j,i) = min(xCenteredSorted);
    end
    
    %x and ylimits
    ymin(j,i) = min([yDataCenteredSorted;yCenteredSorted]);
    
    %x and ylabel
    if j==1
        ylabel('Std of estimates (def)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
    
    %centered x labels
    %set(gca,'ytick',1:4:max(stdPred)+8,'yticklabel',1:4:max(stdPred)+8)
    set(gca,'ytick',1:4:max(max([stdPred;stdData]))+4,...
        'yticklabel',1:4:max(max([stdPred;stdData]))+4)
    set(gca,'xtick',xTickCenteredSorted(3:8:end),...
        'xticklabel',xtickLabel(3:8:end))
    
    %selected legend
    lg = legend(myPlot2(j));
    legend(lg,'boxoff')
end

%x and ylimits
set(h,'ylim',[0 max([stdPred;stdData])])
set(h,'xlim',[min(xmin(:)) - 1  max(xmax(:)) + 1])
SLremoveDeadSpace;
% SLConventionUp(gcf)

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
%get rid of NaN (cmaes doesn't like it)
card = k0(8);
if strcmp(TheModel,'withoutCardinal') || isnan(card)
    k0(8)    = [];
    TheModel = 'withoutCardinal';
end

%fit
%note that function "SLgetLoglCompCmae.m" must be in your code libary
%for cmaes to work
[fitPtmp,negLogl,counteval,stopflag,outputFit] = cmaes('SLgetLoglCompCmae',k0,sqrt(var(k0)),options,data,feature,stimStrength,pstd,priorShape,priorModes,TheModel,varg);

%rearrange
if strcmp(TheModel,'withoutCardinal') || isnan(card)
    fitPdum = nan(10,1);%card
    fitPdum(1:7)  = fitPtmp(1:7);
    fitPdum(9:10) = fitPtmp(8:9);
end
fitPtmp = fitPdum;
