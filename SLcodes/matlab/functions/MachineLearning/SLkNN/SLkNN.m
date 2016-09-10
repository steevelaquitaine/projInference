

% author: Steeve Laquitaine
%   date: 141229
%purpose: Classification with k nearest neighbor method
%
%  usage:
%
%       o = SLkNN(x,y,k,TrainSz)
%       o = SLkNN('Tutorial')
%
%options (varargin)


function o = SLkNN(x,y,k,TrainSz)

%case tutorial
if nargin==1 && strcmp(x,'Tutorial')
    o = SLLoadTextData('Tutorial');
    o.TrainSz = 0.66;
    
    %randomly split data into training and test
    %------------------------------------------
    o.TrainLen = o.TrainSz*o.nlines;
    o.TestLen = o.nlines - o.TrainLen;
    jit = randperm(o.nlines);
    
    %training
    o.quanTrainSet = o.quanData(jit(1:o.TrainLen),:);
    o.qualTrainSet = o.qualData(jit(1:o.TrainLen),:);
    
    %test
    o.quanTestSet = o.quanData(jit(o.TrainLen+1:end),:);
    o.qualTestSet = o.qualData(jit(o.TrainLen+1:end),:);
   
    
    %get neighbors
    %-------------
    
    
end











