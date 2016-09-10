

%author: steeve laquitaine
  %date: 140414
 %usage:e.g., 
     %percentTrain=2/3;
     %i=[{rand(67,288)},{rand(67,288)}]; %2 cells (i.e., classes)
     %paramValsUnq=[45 105]; %One parameter for each class.
     %[datamnTrain,datamnTest,paramValsTrain,paramValsTest]=sortTrainingAndTest(percentTrain,i,paramValsUnq)
%purpose: sort training and test data after concatenating i,a matrix of 
%cells into one matrix.

function [i,iscored,pSettings]=preprocessInstancesAcrossScan(i,pSettings)

%concatenate classes within scan and then concatenate scans.
dscans=[];nScans=length(i);
for iScans=1:nScans
    
    %get scan
    Scan=i{iScans};
    
    %concatenate classes within scan
    d = [];nClasses = length(Scan);
    for iInstance = 1:nClasses
        d(end+1:end+size(Scan{iInstance},1),:) = Scan{iInstance};
        %remember how many instances there originaly was
        nInstances{iScans}(iInstance) = size(Scan{iInstance},1);
    end
    %dscan{iScans}=d;
    
    %now concatenate scans
    dscans(end+1:end+size(d,1),:)=d;
end

%now preprocess
%When z-scoring
%mean across instances is removed (i.e. each voxels response will be 0 
%across all instances)
%std across instances is removed(i.e. each voxels response will have std=1 
%across all instances)
[iscored,pSettings]=preprocessInstances(dscans,pSettings);

%sort back as before in cells of classes within cells of scans
thisRow=1;
i={};
for iScans=1:nScans 
    for iInstance=1:nClasses
        thisEndRow=thisRow+nInstances{iScans}(iInstance) - 1;
        i{iScans}{iInstance}=iscored(thisRow:thisEndRow,:);
        thisRow=thisEndRow+1;
    end
end





