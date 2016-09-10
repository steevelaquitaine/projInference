
%author: steeve laquitaine
  %date: 140411
 %usage: [VarParamSelectivity,Variable]=getParamSelectivity(dataTRAINmn,paramValsTrain)
     %e.g., dataTRAINmn is a matrix (m,n)
     %e.g., paramValsTrain is an arraw of parameters (e.g., directions) 
     %with m values.
     
     %dataTRAINmn=rand(67,288);
     %paramValsTrain=randi(36s0,288,1);
     %[VarParamSelectivity,Variable]=getParamSelectivity(dataTRAINmn,paramValsTrain)
     
     
function [VarParamSelectivity,Variable]=SLgetParamSelectivity(dataTRAINmn,paramValsTrain)

%number of rows (e.g., voxels) and columns (trials)
numrow=size(dataTRAINmn,1);

%Variable
Variable=[1:1:numrow]';

%Parameter values and selectivity
paramValsUnq=unique(paramValsTrain);


%Calculate average data value for each column value (e.g., direction)
%and row's (e.g., voxels) preferred column value, data is averaged over
%repetition of column values. Preferred column value is the column value
%that elicits maximum average data value.
VarParamSelectivity=nan(numrow,1);
parfor row=1:numrow;
    [m,s]=makeStat(dataTRAINmn(row,:)',paramValsTrain);
    [~,pos]=max(m);
    VarParamSelectivity(row)=paramValsUnq(pos);
end
%print message
fprintf('%1s \n','(SLgetParamSelectivity) Getting the selectivity of your variable for a parameter...done')

