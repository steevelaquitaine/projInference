
%SLfMRIforwardModelingAnalysis.m
%
%      Author: Steeve Laquitaine
%        date: 140413 - updated 150108
%     purpose: reconstruct the directions of a motion stimulus from a 
%              pattern of BOLD responses with forward modeling.
%
%       usage: 
%
%             SLfMRIforwardModelingAnalysis('/Users/steeve/data/datafMRI/s02520140220','allMT','myRandomDir','MotionComp','zscore=1','main task','taskNum=2','phaseNum=1')
%
% Description:
%
%     Enter the path of your dataset       
%
%     Train the weights for each voxel and get channels outputs from training
%     BOLD data.
%     
%     Use the weights to get channels output in test data set and compare
%     them to training channels output to reconstruct test displayed directions.
% 
%
%This function requires mrTools (Justin Gardner)

function output = SLfMRIforwardModelingAnalysis(pathDataSet,myROIname,...
    myVariable,curGroup,prepro,scanNumFromDescription,taskNum,phaseNum)

%backup command window
diary Log_SLfMRIforwardModelingAnalysis_140506.txt

%backup m.file with date and time
filename = matlab.desktop.editor.getActiveFilename;
SLBackup(filename)


%BOLD responses "MTBOLDmnAll" (n trials, m voxels) 
%-------------------------------------------------
%Get fMRI and behavioral data from "main task".
cd(pathDataSet)

%Data parameters
%ROI: e.g., Matrix of pooled left and right MT ('allMT',140410) contains 
%67 voxels (col) and 288 trials (8 scans * 36 trials (= 6 dir* 3rep* 2coh)).
%2 trials are discarded at the start of each scan.
%Variable: To get variable names use "getTaskVarnames(v)" after calling 
%viewSet.
%references
%http://gru.brain.riken.jp/doku.php/mrtools/mrtoolssinglepage?s[]=concattseries
v = newView;

%each scan data
v  = viewSet(v,'curGroup',curGroup);
scanNum = viewGet(v,'scanNumFromDescription',scanNumFromDescription,v.curGroup);
i = [];

for j = 1:numel(scanNum)
    
    v = newView;
    v = viewSet(v,'curGroup',curGroup);
    v = viewSet(v,'curScan',scanNum(j));
        
    %detrend and high-pass filter the scans you select, get rid of any 
    %junk frames, do a de-trending step using a hi-pass filter if you 
    %specify that, and then concatenate the resulting time series into a 
    %single time series that will be saved into a new group called 
    %Concatenation. The mat files are saved in folder "concatenation".
    [v,params] = concatTSeries(v,[],'justGetParams=1','defaultParams=1',...
        'scanList',scanNum(j));
    v = concatTSeries(v,params);
    fprintf('\n %12s \n',['Scan ',num2str(scanNum(j)),'Concatenation...done'])
    
    %update group and scan to concatenated data
    v = newView;
    v = viewSet(v,'curGroup',params.newGroupName);
    nScans = viewGet(v,'nScans');
    v = viewSet(v,'curScan',nScans);

    %extract data
    myROI = loadROITSeries(v,myROIname);
    myROI =getInstanceMatrix(v,myVariable,myROI,taskNum,phaseNum,'n=inf');
    myROI = myROI{1};
    output.i{j} = myROI.classify.instances;
    output.paramValsUnq(:,j) = unique(myROI.classify.paramVals);
    output.paramsVal{j} = myROI.classify.paramVals;
end
fprintf('\n %12s \n','Get BOLD responses...done')


%----------------------------------------------------------
%%Preprocessing: Z-score "MTBOLDmnAll" (n trials, m voxels) 
%-----------------------------------------------------------
%Matrix must be one row for each instance, one column for each voxel
%so d=nxk where n is number of instances and m is number of voxels
%(instances (i.e., voxels responses at each trial) are z-scored across 
%trials for each voxel independently of the classes (i.e.,coherence or 
%motion direction).i.e., mean and std voxels responses ?calculated for 
%each voxel (vector of m means) and substracted (/divided) to each 
%instance of voxel response.
%If z-scoring works fine std and mean(iMatrix)=1 and 0 respectively.
[output.iscored,iMatrix,pSettings] = preprocessInstancesAcrossScan(output.i,prepro);
fprintf('\n %12s \n','Z-score of BOLD responses...done')


%-------------------------------------------------------------------
%%(Optional) Decode parameters in training dataset to check that 
%weights are good.
%-------------------------------------------------------------------
% %parameters
% %put everyting in training 
% [BOLDmnTrain,~,paramValsTrain]=sortTrainingAndTestPerScan('numScanTrain=8',...
%     output.iscored,...
%     output.paramsVal);
% 
% %get voxels' selectivity for motion direction
% [VarParamSelectivity,Var]=SLgetParamSelectivity(BOLDmnTrain,paramValsTrain);
% 
% %sort voxels by selectivity and draw. Voxels are sorted in the training 
% %data and this sorting is applied to the test data set. Code was checked
% %with duymmy data and works.
% close all
% [BOLDmnTrainSorted,paramSelectivity,VarSorted]=SLsortByParamSelectivity(BOLDmnTrain,...
%     paramValsTrain,VarParamSelectivity);
% 
% 
% % FE and decode training parameters
% [decodedParamTrain,~,Rdec,Pdec]=SLdoFEandReconstruct(BOLDmnTrain,BOLDmnTrain,...
%     paramValsTrain,paramValsTrain,paramSelectivity,'display=on');


%-----------------------------------
%Reconstruct classes in test dataset 
%-----------------------------------
%Use cross-validation: leave-one-scan out test on remaining scan.
%Make training and Test dataset. Training dataset is 2/3 of the 
%dataset for each class. e.g. if a class (e.g., direction 45 degrees) was 
%repeated 48 times we take 32 repetitions of direction 45 out of 48. Same
%for all directions.

%scans
output.numScans = length(output.iscored);
output.numScanTrain = 7;
output.numScansTest = output.numScans - output.numScanTrain;

%Order scans for training '1 : numScanTrain' and test 'end-numScanTrain'
output.TrainingSetList = nchoosek(1:8,output.numScanTrain);

%sets for training and decoding
output.numSets = size(output.TrainingSetList,1);
TestSetList = nan(output.numSets,output.numScansTest);
for j = 1 : output.numSets
    TestSetList(j,:) = setdiff(1 : output.numScans,output.TrainingSetList(j,:));
end
SetsList = [output.TrainingSetList TestSetList];

%decode parameters with each set
reconstParamtestBkup=[];
paramValsTestBkup=[];
RrecBkup = [];
PrecBkup = [];

for j = 1 : output.numSets
    
    %update
    fprintf('\n %12s \n',['Decoding test set ',num2str(j)])
    
    %scans data ordered to be sorted in training and test data
    ThisSetList = SetsList(j,:);    
    [BOLDmnTrain,BOLDmnTest,paramValsTrain,paramValsTest] = sortTrainingAndTestPerScan(['numScanTrain=',num2str(output.numScanTrain)],...
        output.iscored(ThisSetList),...
        output.paramsVal(ThisSetList));
    
    %get voxels' selectivity for motion direction
    [VarParamSelectivity,Var] = SLgetParamSelectivity(BOLDmnTrain,paramValsTrain);
    
    %sort voxels by selectivity and draw. Voxels are sorted in the training
    %data and this sorting is applied to the test data set. Code was checked
    %with duymmy data and works.
    close all
    [BOLDmnTrainSorted,paramSelectivity,VarSorted] = SLsortByParamSelectivity(BOLDmnTrain,...
        paramValsTrain,VarParamSelectivity,'display=off');
    xlabel('Motion directions (degrees)');
    ylabel('Voxels by p referred direction (degrees)')
        
    %apply the same sorting to test data set (if there are test data)
    if isempty(BOLDmnTest)==1
        BOLDmnTest = BOLDmnTrain;
        paramValsTest = paramValsTrain;
        fprintf('\n %12s \n','Same data (all) are used for training and testing')
    end     
    BOLDmnTestSorted=BOLDmnTest(VarSorted,:);
    
    %FE and reconstruction of test parameters
    reconstParamtest = SLdoFEandReconstruct(BOLDmnTrainSorted,BOLDmnTestSorted,...
        paramValsTrain,paramValsTest,paramSelectivity,'display=off');
    
    %keep track of test and reconstructed parameters.
    %parameters
    reconstParamtestBkup(:,j) = reconstParamtest;
    paramValsTestBkup(:,j) = paramValsTest;
end
fprintf('\n %12s \n','Cross-validated reconstruction done')


%calculate decoding accuracy
%get all parameters and reconstructed parameter
reconstParamtestAll=reconstParamtestBkup(:);
paramValsTestAll=paramValsTestBkup(:);
paramValsUnq=unique(paramValsTestAll);
AccuracyThisParamVal=nan(length(paramValsUnq),1);
numberOfMatch=nan(length(paramValsUnq),1);
numberOfRepParamsValue=nan(length(paramValsUnq),1);

for j = 1 : length(paramValsUnq)
    
    %get each test parameter value
    thisParamVal = paramValsUnq(j);
    
    %get reconstructed parameters for this test parameter
    reconstParamForThisParamVal = reconstParamtestAll(paramValsTestAll==thisParamVal);
    numberOfRepParamsValue(j) = numel(reconstParamForThisParamVal);
    
    %calculate probability that reconstructed value matches test parameter
    numberOfMatch(j) = sum(reconstParamForThisParamVal==thisParamVal);
    AccuracyThisParamVal(j) = numberOfMatch(j)/numberOfRepParamsValue(j);
end

%overall decoding accuracy (probability and correlation)
Accuracy=sum(numberOfMatch)/sum(numberOfRepParamsValue);
[Rrec,Prec]=corrcoef(reconstParamtestAll,paramValsTestAll);

%Plot mean reconstructed directions with decoding information
%number of voxels
numVoxels=size(BOLDmnTestSorted,1);

%information about the task
%number of instances for each classes
numInstPerClass=size(output.iscored{1}{1},1);

%number of classes
numClasses=size(output.iscored{1},2);


%------------------------
%draw reconstructed class
%------------------------
hold all
[means,stds] = drawCircStat(reconstParamtestAll,paramValsTestAll,[],[],'errorbar','clean=off');
plot([1 360],[225 225],'b--')
title({['Cross-validated decoding of motion directions',date],...
    ['Accuracy ',num2str(Accuracy),...
    ', R=',num2str(Rrec(2,1)),', p=',num2str(Prec(2,1)),'(Pearson)'],...
    ['(',num2str(output.numScanTrain),' train/',...
    num2str(output.numScansTest),' test scans,',...
    num2str(numInstPerClass),' trials/class, ',...
    num2str(numClasses)',' classes, ',...
    'ROI: ',myROIname,', ',num2str(numVoxels),' voxels, ',...
    'Preprocess with ',prepro,'pca:',num2str(pSettings.pca),')']},...
    'fontsize',11)
ylabel('Mean decoded directions (degrees)','fontsize',14)
xlabel('Motion directions in Test (degrees)','fontsize',14)


% %% draw accuracies for different ROIs
% accuracies=[0.30903 0.21];
% ROIs={'V1toMT','MT'};
% figure;
% hold all
% SLdrawBar(accuracies,ROIs,[1 2])
% plot([1:numel(accuracies)],1./[numClasses numClasses],'r--')
% ylabel('Accuracies','fontsize',14)
% title(date)


% %% draw BOLD responses for each motion direction sorted by voxels' 
% %preference for direction
% %get BOLD responses (m voxels, n instances)
% output.numScans=length(output.iscored);
% [BOLDmnAll,~,paramValsAll]=sortTrainingAndTestPerScan(['numScanTrain=',num2str(output.numScans)],...
%     output.iscored,...
%     output.paramsVal);
% fprintf('\n %12s \n','Get all BOLD responses...done')
% 
% %sort voxels by selectivity
% %get selectivity
% [VarParamSelectivity,Var]=SLgetParamSelectivity(BOLDmnAll,paramValsAll);
% fprintf('\n %12s \n','Get voxels selectivity...done')
% 
% %sort voxels by selectivity and draw. Voxels are sorted in the training
% %data and this sorting is applied to the test data set. Code was checked
% %with duymmy data and works.
% [BOLDmnSorted,paramSelectivity,VarSorted]=SLsortByParamSelectivity(BOLDmnAll,...
%     paramValsAll,VarParamSelectivity);
% xlabel('Motion directions (degrees)');
% ylabel('Voxels by p referred direction (degrees)')
% fprintf('\n %12s \n','Sort voxels by selectivity...done')
% 
% %plot mean BOLD response (across repetitions) for each motion direction
% %with necessary informations
% paramsValuniq=unique(paramValsAll);
% numparamsValuniq=numel(paramsValuniq);
% mall=nan(size(BOLDmnSorted,1),numparamsValuniq);
% sall=nan(size(BOLDmnSorted,1),numparamsValuniq);
% for j=1:size(BOLDmnSorted,1)
%     [m,s]=makeStat(BOLDmnSorted(j,:)',paramValsAll);
%     mall(j,:)=m;
%     sall(j,:)=s;
% end
% figure ('color','w')
% imagesc(mall)
% colorbar
% 
% %number of classes
% numClasses=size(output.iscored{1},2);
% 
% %number of voxels
% numVoxels=size(BOLDmnSorted,1);
% 
% set(gca,'ytick',1:60:numel(paramSelectivity),'yticklabel',paramSelectivity(1:60:end))
% set(gca,'xtick',1:numparamsValuniq,'xticklabel',paramsValuniq)
% ylabel({'Voxels sorted by','direction selectivity (degrees)'},'fontsize',14)
% title({['Bold responses (m voxels, motion directions)',date],...
%     ['(',num2str(output.numScans),' scans, ',...
%     num2str(numClasses)',' classes, ',...
%     'ROI: ',myROIname,', ',num2str(numVoxels),' voxels, ',...
%     'Preprocess with ',prepro,'pca:',num2str(pSettings.pca),')']},...
%     'fontsize',11)
% colormap('hot')
% axis square
% xlabel('Motion directions (degrees)','fontsize',14)
% fprintf('\n %12s \n','Plot BOLD responses with voxels sorted by selectivity...done')


%backup and end diary
save([filename(1:end-2),'.mat'])
diary off




