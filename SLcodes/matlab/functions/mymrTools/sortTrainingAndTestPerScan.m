
%author: steeve laquitaine
  %date: 140411
 %usage:e.g., 
     %percentTrain=2/3;
     %i=[{rand(67,288)},{rand(67,288)}]; %2 cells (i.e., classes)
     %paramValsUnq=[45 105]; %One parameter for each class.
     %[datamnTrain,datamnTest,paramValsTrain,paramValsTest]=sortTrainingAndTestPerScan(percentTrain,i,paramValsUnq)
%purpose: sort training and test data after concatenating i,a matrix of 
%cells into one matrix.
     
function [datamnTrain,datamnTest,paramValsTrain,paramValsTest]=sortTrainingAndTestPerScan(numScanTrain,i,paramVals)

%parameters
%number of scans
numScans=length(i);

%number of training scans
numScanTrain=str2double(numScanTrain(find(numScanTrain=='=')+1:end));

%number of test scans
numScanTest=numScans-numScanTrain;

%%TRAINING
%concatenate classes within scan and then concatenate scans.
%e.g., 8 scans x 6 classes x 6 repetitions becomes a matrix of 288 trials 
%x 67 voxels.
paramValsTrain=[];
dscans=[];
TrainScanNum=[];
for iScans=1:numScanTrain
    
    %set train scan
    TrainScanNum(iScans)=iScans;
    
    %get scan
    ScanTrain=i{TrainScanNum(iScans)};
    
    %concatenate classes within scan
    d = [];nClasses = length(ScanTrain);
    for iInstance = 1:nClasses
        d(end+1:end+size(ScanTrain{iInstance},1),:) = ScanTrain{iInstance};
        %remember how many instances there originaly was
        nInstancesTrain{TrainScanNum(iScans)}(iInstance) = size(ScanTrain{iInstance},1);
    end
    
    %now concatenate scans
    dscans(end+1:end+size(d,1),:)=d;
    
    %concatenate parameter values
    paramValsTrain(end+1:end+size(d,1),:)=paramVals{TrainScanNum(iScans)};
end

%matrices of data (m voxels, n-ntraining trials) 
datamnTrain=dscans';


%%TEST
%concatenate remaining scans
dscans=[];
TestScanNum=[];
paramValsTest=[];

%if there is a test scan
if numScanTest>0
    for iScans=1:numScanTest
        
        %set test scan
        TestScanNum(iScans)=numScanTrain+iScans;
        
        %get scan
        ScanTest=[];
        ScanTest=i{TestScanNum(iScans)};
        
        %concatenate classes within scan
        d=[];nClasses=length(ScanTest);
        for iInstance = 1:nClasses
            d(end+1:end+size(ScanTest{iInstance},1),:) = ScanTest{iInstance};
            %remember how many instances there originaly was
            nInstancesTest{iScans}(iInstance) = size(ScanTest{iInstance},1);
        end
        
        %now concatenate scans
        dscans(end+1:end+size(d,1),:)=d;
                
        %concatenate parameter values
        paramValsTest(end+1:end+size(d,1),:)=paramVals{TestScanNum(iScans)};
    end
end

%matrices of data (m voxels, end-n trials) 
datamnTest=dscans';




% %concatenate training scan
% for j=1:numScanTrain
%     
%     %get data and parameter values for each scan
%     %training 
%     numI=size(iScan{j},1);
%     numItrain=round(percentTrain*numI);
%     itrain{j}=i{j}(1:numItrain,:);
%     itrainParamVal{j}=repmat(paramValsUnq(j),numItrain,1);
%     
%     %test
%     itest{j}=i{j}(numItrain+1:numI,:);
%     itestParamVal{j}=repmat(paramValsUnq(j),numI-numItrain,1);
% 
%     %concatenate in matrices (n rows,m colfire) and arrays of
%     %parameter values
%     itrainAll=[itrainAll; itrain{j}];
%     itestAll=[itestAll; itest{j}];
%     paramValsTrain=[paramValsTrain; itrainParamVal{j}];
%     paramValsTest=[paramValsTest; itestParamVal{j}];
% end
% 
