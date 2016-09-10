

%slfMRIcheckVMfitTuningSignif.m
%
% author: steeve laquitaine
%purpose: check significance of observed voxel vm fit tuning strength 
%         against permuted null tunings calculated from fitting von 
%         mises to permutation of voxel responses.
%
%  usage: 
%
%       [sig,nVecs,actualnVec,actualnVecStrng,meantunStrngBydir,tunStrngCIbyDir,vecdir] = slfMRIcheckVMfitTuningSignif(y,x,nperm)
%
% 

function [sig,nVecs,actualnVec,actualnVecStrng,meantunStrngBydir,tunStrngCIbyDir,vecdir,minsqe] = slfMRIcheckVMfitTuningSignif(y,x,nperm)

%get voxel tuning directions and strength by permutation
[nVecs,actualnVec,minsqe] = slfmriGetSingleVoxNullvmfitTuningDist(y,x,nperm);

%get if actual data vector strength at data tuned direction is outside the
%confidence interval of the strength of the same direction predicted by 
%the permuted data (see if significanty sharper than predicted by chance)
%calculate mean and CI of tuning strength by tuning direction 
[meantunStrngBydir,vecdir] = slCalculateMeanSortedBy1Var(nVecs.deg.tunStrngth',nVecs.deg.mean);
tunStrngCIbyDir = slCalculateCISortedBy1Var(nVecs.deg.tunStrngth',nVecs.deg.mean);
actualnVec.deg.mean = actualnVec.datadegmean;
ActtunDir_loc = vecdir==round(actualnVec.deg.mean);

%if actual tuning not produced by permutations, 
%pick the closest value
if sum(ActtunDir_loc)==0
    [~,ActtunDir_loc] = min(abs(vecdir - round(actualnVec.deg.mean)));    
end
sig = actualnVec.tunStrngth > tunStrngCIbyDir(ActtunDir_loc,2);

%backup
actualnVecStrng = actualnVec.tunStrngth;
