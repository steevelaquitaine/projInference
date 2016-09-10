
%author: steeve laquitaine
  %date: 140411
 %usage:e.g., 
     %percentTrain=2/3;
     %i=[{rand(67,288)},{rand(67,288)}]; %2 cells (i.e., classes)
     %paramValsUnq=[45 105]; %One parameter for each class.
     %[datamnTrain,datamnTest,paramValsTrain,paramValsTest]=sortTrainingAndTest(percentTrain,i,paramValsUnq)
%purpose: sort training and test data after concatenating i,a matrix of 
%cells into one matrix.
     
function [datamnTrain,datamnTest,paramValsTrain,paramValsTest]=sortTrainingAndTest(percentTrain,i,paramValsUnq)

itrainAll=[];
itestAll=[];
paramValsTrain=[];
paramValsTest=[];
for j=1:length(i)
    
    %get data and parameter values for each class
    %training 
    numI=size(i{j},1);
    numItrain=round(percentTrain*numI);
    itrain{j}=i{j}(1:numItrain,:);
    itrainParamVal{j}=repmat(paramValsUnq(j),numItrain,1);
    
    %test
    itest{j}=i{j}(numItrain+1:numI,:);
    itestParamVal{j}=repmat(paramValsUnq(j),numI-numItrain,1);

    %concatenate in matrices (n rows,m colfire) and arrays of
    %parameter values
    itrainAll=[itrainAll; itrain{j}];
    itestAll=[itestAll; itest{j}];
    paramValsTrain=[paramValsTrain; itrainParamVal{j}];
    paramValsTest=[paramValsTest; itestParamVal{j}];
end

%matrices of data (m voxels,n trials) for forward encoding model
datamnTrain=itrainAll';
datamnTest=itestAll';