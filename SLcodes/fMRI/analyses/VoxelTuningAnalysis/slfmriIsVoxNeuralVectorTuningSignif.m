

%author: steeve laquitaine
%purpose : check whether a voxel is significantly tuned for a particular 
%          motion direction : whether neural vector is significantly different
%          from the null distribution of neural vectors obtained from permutating
%          voxel responses
%
%   usage: 
%           nVec = slfmriGetSingleVoxNullneuralVectorDist([10 1 4 7],[0 90 180 270],1000)
%           [o,actualnVecLen,ci] = slfmriIsVoxNeuralVectorTuningSignif(nVecs,actualnVec);



function [o,actualnVecLen,ci] = slfmriIsVoxNeuralVectorTuningSignif(nVecs,actualnVec)

%get voxel null tuning vector length confidence interval
parfor i = 1 : length(nVecs)       
    nVecLen(i) = slGetVectorLength(nVecs(i).coord.mean(1),nVecs(i).coord.mean(2));
end

%get confidence interval and 
%take care of numerical imprecision
[~,ci] = slMakeCI(nVecLen,.95);
ci = round(ci*100)/100;

%get voxel actual tuning vector length
actualnVecLen = slGetVectorLength(actualnVec.coord.mean(1),actualnVec.coord.mean(2));
actualnVecLen = round(actualnVecLen*100)/100;

%check whether the actual vector length is significant
%(outside the confidence interval)
o = slIsBetween(actualnVecLen,ci(1),ci(2));
o = 1 - o;