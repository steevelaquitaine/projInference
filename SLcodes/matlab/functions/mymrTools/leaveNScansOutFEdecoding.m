
     %Author: Steeve Laquitaine
       %date: 140421
  %file name: leaveNScansOutFEdecoding.m
      %usage: 
      %e.g., 
      %iscored=repmat({repmat({randi(10,6,67)},1,6)},1,8);
      %stimVal=repmat({randi(360,36,1)},1,8);
      %[Accuracy,Rdec,Pdec]=leaveNScansOutFEdecoding(iscored,7,stimVal);
      
      %8 scans, 6 stimuli classes (e.g., 6 motion directions), 6
      %stimuli repetitions/scan (36 trials), 67 voxels in the ROI, 7 scans
      %for training, stimVal are presented stimuli for each instance.
   
    %purpose: decode by training forward modeling on N training scans and
              %decode stimulus in end-N remaining scans.

%description:
    %1. Train the forward modeling weights for each voxel and get channels 
    %outputs from training BOLD data.
    %2 use the weights to get channels output in test data set and compare
    %them to training channels output to decode test displayed stimulus.
    
    %Use cross-validation: leave-one-scan out test on remaining scan.
    %Make training and Test dataset. Training dataset is 2/3 of the 
    %dataset for each class. e.g. if a class (e.g., direction 45 degrees) was 
    %repeated 48 times we take 32 repetitions of direction 45 out of 48. Same
    %for all directions.

function [Accuracy,Rdec,Pdec]=leaveNScansOutFEdecoding(iscored,numScanTrain,stimsVal,myROIname,pSettings,OptionDisp)

%number of scans
numScans=length(iscored);

%number of test scans
numScansTest=numScans-numScanTrain;

%order scans for training '1:numScanTrain' and test 'end-numScanTrain'
TrainingSetList=nchoosek(1:8,numScanTrain);

%create sets for training and decoding
numSets=size(TrainingSetList,1);
TestSetList=nan(numSets,numScansTest);
for j=1:numSets
    TestSetList(j,:)=setdiff(1:numScans,TrainingSetList(j,:));
end
SetsList=[TrainingSetList TestSetList];

%decode stimuli via cross-validation
decodedStimtestBkup=[];
stimValTestBkup=[];
RrecBkup=[];
PrecBkup=[];

for j=1:numSets
    
    %update
    fprintf('\n %12s \n',['(leaveNScansOutFEdecoding) Decoding test set ',num2str(j)])
    
    %get scans data ordered to be sorted in training and test data
    ThisSetList=SetsList(j,:);    
    [BOLDmnTrain,BOLDmnTest,stimValTrain,stimValTest]=sortTrainingAndTestPerScan(['numScanTrain=',num2str(numScanTrain)],...
        iscored(ThisSetList),...
        stimsVal(ThisSetList));
    
    %get voxels' selectivity for motion direction
    [VarStimSelectivity,~]=SLgetParamSelectivity(BOLDmnTrain,stimValTrain);
    
    %sort voxels by selectivity and draw. Voxels are sorted in the training
    %data and this sorting is applied to the test data set. Code was checked
    %with duymmy data and works.
    %close all
    [BOLDmnTrainSorted,stimSelectivity,VarSorted]=SLsortByParamSelectivity(BOLDmnTrain,...
        stimValTrain,VarStimSelectivity,'OptionDisp=0');
    xlabel('Motion directions (degrees)');
    ylabel('Voxels by p referred direction (degrees)')
        
    %apply the same sorting to test data set (if there are test data)
    if isempty(BOLDmnTest)==1
        BOLDmnTest=BOLDmnTrain;
        stimValTest=stimValTrain;
        fprintf('\n %12s \n','Same data (all) are used for training and testing')
    end     
    BOLDmnTestSorted=BOLDmnTest(VarSorted,:);
    
    %FE and decoding of test stimuli
    decodedStimtest=SLdoFEandReconstruct(BOLDmnTrainSorted,BOLDmnTestSorted,...
        stimValTrain,stimValTest,stimSelectivity,'display=off');
    
    %store test and decoded stimuli
    decodedStimtestBkup(:,j)=decodedStimtest;
    stimValTestBkup(:,j)=stimValTest;
end
fprintf('\n %12s \n','(leaveNScansOutFEdecoding) Cross-validated decoding of test data set....done')


%calculate decoding accuracy
%get all stimulus and decoded parameter
decodedStimtestAll=decodedStimtestBkup(:);
stimValTestAll=stimValTestBkup(:);

stimValUnq=unique(stimValTestAll);
AccuracyThisStimVal=nan(length(stimValUnq),1);
numberOfMatch=nan(length(stimValUnq),1);
numberOfRepStimsValue=nan(length(stimValUnq),1);
for j=1:length(stimValUnq)
    
    %get each test parameter value
    thisStimVal=stimValUnq(j);
    
    %get decoded stimulus for this test parameter
    decodedStimForThisStimVal=decodedStimtestAll(stimValTestAll==thisStimVal);
    numberOfRepStimsValue(j)=numel(decodedStimForThisStimVal);
    
    %calculate probability that decoded value matches test parameter
    numberOfMatch(j)=sum(decodedStimForThisStimVal==thisStimVal);
    AccuracyThisStimVal(j)=numberOfMatch(j)/numberOfRepStimsValue(j);
end

%overall decoding accuracy (probability that decoded and actual match). 
Accuracy=sum(numberOfMatch)/sum(numberOfRepStimsValue);

%use circular correlation for circular data.
Rdec=nan(2,2);
Pdec=nan(2,2);
%[Rdec,Pdec]=corrcoef(decodedStimtestAll,stimValTestAll);
Rdec=Rdec(2,1);
Pdec=Pdec(2,1);

%Plot mean decoded directions with decoding information
%number of voxels
numVoxels=size(BOLDmnTestSorted,1);

%information about the task
%number of instances for each class
numInstPerClass=size(iscored{1}{1},1);

%number of classes of stimuli
numClasses=size(iscored{1},2);

%draw decoded stimuli against actual stimuli
if strcmp(OptionDisp,'display=on')
    hold all
    drawCircStat(decodedStimtestAll,stimValTestAll,[],[],'errorbar');
    plot([1 360],[225 225],'b--')
    title({['(leaveNScansOutFEdecoding) Cross-validated decoding of motion directions',date],...
        ['Accuracy ',num2str(Accuracy),...
        ', R=',num2str(Rdec),', p=',num2str(Pdec),'(Pearson)'],...
        ['(',num2str(numScanTrain),' train/',...
        num2str(numScansTest),' test scans,',...
        num2str(numInstPerClass),' trials/class, ',...
        num2str(numClasses)',' classes, ',...
        'ROI: ',myROIname,', ',num2str(numVoxels),' voxels, ',...
        'Preprocess with zscore:',num2str(pSettings.zscore),' pca:',num2str(pSettings.pca),')']},...
        'fontsize',11)
    ylabel('Mean decoded directions (degrees)','fontsize',14)
    xlabel('Motion directions in Test (degrees)','fontsize',14)
end
