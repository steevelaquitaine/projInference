

%SLfitLearningBayesianModel.m
%
%
%     author: steeve laquitaine
%       date: 160608
%     status: Complete
%    purpose: Model motion direction estimation data as Bayesian
%             Also fit alternative (suboptimal Bayesian models)
%
%
%      usage:
%
%           SLfitLearningBayesianModel({'sub01'},[23.78 7.56 1.88 0.66 1.71 2.77 5.30 NaN 0.001 83.63 NaN 20],'experiment','vonMisesPrior','filename','datafit','MAPReadout','MaxLikelihoodFit','fminsearch');
%
%
%inputs:
%
%         initp: 3 kl, 4 kp, 1 kcardinal, 1 prand, 1 Kmotor, 1 tail coefficient, 1 learning tau
%
%nested functions : slfitMaxLLH

function [fitP,fitPbkp,R2,sdata,fitPt,negLogl,negLoglbkp,...
    Logl_pertrialBestfit,output] = SLfitLearningBayesianModel(subjects,...
    initp,...
    varargin)

%time
t0 = tic;

%help
if ieNotDefined('subjects')
    help SLfitLearningBayesianModel
    return
end

%initialize outputs
output = [];
output.code = 'SLfitLearningBayesianModel.m';

%data path
vrg = varargin;
if slIsInput(vrg,'vonMisesPrior')
    
    %get data
    if slIsInput(vrg,'dataPathVM')
        dataPathVM = varargin{find(strcmp(varargin,'dataPathVM'))+1};
        cd(dataPathVM)
    else
        fprintf('%s \n','(slfitLearningBayes) You need to set dataPath ...')
        dataPathVM ='~/data/dataPsychophy/proj01_priorStrength';%uigetdir(cd,'Pickup your project e.g., /dataPsychophy/Exp01...');
        cd(dataPathVM)
    end
    vararginVM = [varargin,'dataPath',dataPathVM];
    
    %Load or create database
    if sum(strcmp(varargin,'LoadDataBase'))==0
        fprintf('%s \n','(slfitLearningBayes) Creating database ...')
        databankVM = SLMakedatabank(subjects,vararginVM);
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

%backup file name
if sum(strcmp(varargin,'filename'))
    filename = varargin{find(strcmp(varargin,'filename'))+1};
else
    filename  = input('SLfitLearningBayesianModel) Please set the name of the .mat file that will be saved. (e.g., dataSub01)','s');
end

%check initp
if size(initp) < 12
    fprintf('%s \n','(SLfitLearningBayesianModel) Not enough initial parameters ...')
    keyboard
end

%(case cardinal prior or not)
if isnan(initp(8)); TheModel='withoutCardinal'; else TheModel='withCardinal'; end

%Initialize outputs
R2 =[]; udata =[]; sdata =[]; dataOut =[]; predOut =[]; FAbp=[];
FA=[];fitPt =[]; fitPbkp =[]; negLoglbkp=[]; Logl_pertrialBestfit=[];
negLogl=[];

%You can set fit parameters (or not,'[]').
fitP = [];

%check analysis
posTarget = SLfindword(varargin,{'MaxLikelihoodFit'});
if sum(posTarget)==0
    fprintf('%s \n',['(SLfitBayesianModel) You need to input an analysis:',...
        ' Possible analyses are: '],...
        ' - MaxLikelihoodFit, ')
end

%Maximum likelihood fit
if sum(strcmp(varargin,'MaxLikelihoodFit'))
    fprintf('%s \n','(slfitMaxLLH) Now fitting the model with max LLH')
    %fit
    
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
        slfitMaxLLH(databankVM,fitP,initp,filename,output,varargin);
    
    %save data
    output.duration = toc(t0);
    mkdir(['fitLearningBayes' subjects{1}])
    cd(['fitLearningBayes' subjects{1}])
    slPrintfStr('slfitLearningBayesianModel',['Saving data in' pwd])
    save(filename)
end

%(If we want data mean, std and distribution with models' predictions)
if sum(strcmp(varargin,'modelPredictions'))==1
    
    %von Mises prior
    if sum(strcmp(varargin,'vonMisesPrior'))
        databank = databankVM;
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
        
        %# of subjects
        inputSubjects = subjects;
        subjectsNum = unique(databank.subjects);
        numSub = numel(subjectsNum);
        
        %data and predictions will be averaged across subjects
        if numSub > 1
            fprintf(['(SLdrawModelsPredictionCentered) More that 1 subject have been found',...
                '. Data and models predictions will be averaged across subjects /n'])
        end
        
        %all experimental conditions
        %case von Mises prior
        %--------------------
        output.uniqCond = SLuniqpair([databank.Pstd databank.stimStrength ...
            databank.stimFeatureDeg]);
        numcond = size(output.uniqCond,1);
        
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
        
        %subjects
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
            [meanDatatmp,stdDatatmp,dataCondtmp] = makeDataMeanAndStd(data,...
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
                
                %store subjects' best fit parameters
                output.fitP(sub,:) = fitP.p;
                
                %add NaN weightFatTail parameter to old data
                if length(output.fitP(sub,:)) < 11
                    output.fitP(sub,end+1) = NaN;
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
        end
        output.meanDataOvSub = SLmakeColumn(output.meanDataOvSub);
        output.stdDataOvSub  = nanmean(output.stdData,2);
        output.meanDisDataOvSub = nanmean(output.PdisData,3);
        
        %calculate predictions mean and std and distributions
        for i = 1 : size(output.meanPred,1)
            PredCircMeanOvSub = SLcircMeanStd(output.meanPred(i,:)','polar');
            output.meanPredOvSub(i) = PredCircMeanOvSub.deg.mean;
        end
        output.meanPredOvSub = SLmakeColumn(output.meanPredOvSub);
        output.stdPredOvSub  = nanmean(output.stdPred,2);
        output.meanDisPredOvSub = nanmean(output.PdisPred,3);
        
        %draw data mean and std and their models' predictions
        SLdrawModelsPredictionCentered(output.meanDataOvSub,output.stdDataOvSub,...
            output.meanPredOvSub,output.stdPredOvSub,...
            output.uniqCond,priorModes,priorShape,'yCentered')
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
    fprintf('%s \n',['(slfitMaxLLH) Now fitting the model with',...
        'maximum likelihood method to get best fit parameters standard error'])
    %fit
    [negLogl,negLoglbkp,fitP,fitPbkp,Logl_pertrialBestfit,output] = ...
        slfitMaxLLH(databankVM,fitP,initp,filename,output,varargin);
    %save
    save(filename)
end

%time
output.fitDuration = toc(t0);

%Make sure everything's saved in this directory
if length(subjects)==1
    mkdir(['slfitLearningBayes/' subjects{1}])
    cd(['slfitLearningBayes/' subjects{1}])
    exp = varargin{find(strcmp(varargin,'experiment'))+1};
    save(['datafit', subjects{1}(end-1:end),'_',exp,'_slfitLearningBayes'],'output')
    fprintf('%s %s %s \n','(slfitLearningBayesianModel) Saving fit results in', ['"data/slfitLearningBayes/',subjects{1}],'"')
end



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
