

%slfMRIplotVoxTuningVsChance.m
%
% author: steeve laquitaine
%purpose: check significance of observed tuning vector against permuted 
%         null vectors calculated for a tuning function with 
%         Georgopoulos population vector analysis (by shuffling responses 
%         and associated variable)
%
%  usage: 
%
%       sig = slfMRIcheckTuningSignif(y,x,nperm)
%
% 


function [sig,nVecs,actualnVec,actualnVeclen,meanVlenBydir,VlenCIbyDir,vecdir] = slfMRIcheckTuningSignif(y,x,nperm)

%get neural vector directions and lengths by permutation
[nVecs,actualnVec] = slfmriGetSingleVoxNullneuralVectorDist(y,x,nperm);
nVecsXmean = nVecs.coord.mean(:,1);
nVecsYmean = nVecs.coord.mean(:,2);
nullVeclen = nan(length(nVecs),1);
parfor i = 1 : length(nVecs.deg.mean)    
    nullVeclen(i) = slGetVectorLength(nVecsXmean(i),nVecsYmean(i));   
end
actualnVeclen = slGetVectorLength(actualnVec.coord.mean(1),actualnVec.coord.mean(2));

% get if actual data vector length at data vector direction is outside the
%confidence interval of the length of the same direction predicted by 
%the permuted data (see if significanty sharper than predicted by chance)
%calculate mean and CI of vector length by vector direction 
[meanVlenBydir,vecdir] = slCalculateMeanSortedBy1Var(nullVeclen',nVecs.deg.mean);
VlenCIbyDir = slCalculateCISortedBy1Var(nullVeclen',nVecs.deg.mean);
ActVDir_loc = vecdir==round(actualnVec.deg.mean);

%if actual vector not produced by permutations, 
%pick the closest produced value
if sum(ActVDir_loc)==0
    [substitutenvecdir,ActVDir_loc] = min(abs(vecdir - round(actualnVec.deg.mean)));    
end

sig = actualnVeclen < VlenCIbyDir(ActVDir_loc,1) | actualnVeclen > VlenCIbyDir(ActVDir_loc,2)
