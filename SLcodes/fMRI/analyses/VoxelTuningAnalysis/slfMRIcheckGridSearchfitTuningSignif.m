

%slfMRIcheckGridSearchfitTuningSignif.m
%
% author: steeve laquitaine
%purpose: check significance of observed voxel vm fit tuning strength 
%         against permuted null tunings calculated from fitting von 
%         mises to permutation of voxel responses.
%
%  usage: 
%
%       [sig,nVecs,actualnVec,actualnVecStrng,meantunStrngBydir,tunStrngCIbyDir,vecdir] = slfMRIcheckGridSearchfitTuningSignif(y,x,nperm)
%
% 

function [sig,nVecs,actualnVec,actualnVecStrng,meantunStrngBydir,tunStrngCIbyDir,vecdir,minsqe] = slfMRIcheckGridSearchfitTuningSignif(y,x,nperm)

%get voxel tuning directions and strength by permutation
[nVecs,actualnVec,minsqe] = slfmriGetSingleVoxNullgridSearchfitTuningDist(y,x,nperm);

%get if actual data vector strength at data tuned direction is outside the
%confidence interval of the strength of the same direction predicted by 
%the permuted data (see if significanty sharper than predicted by chance)
%calculate mean and CI of tuning strength by tuning direction 
[meantunStrngBydir,vecdir] = slCalculateMeanSortedBy1Var(nVecs.deg.tunStrngth',nVecs.deg.mean);
tunStrngCIbyDir = slCalculateCISortedBy1Var(nVecs.deg.tunStrngth',nVecs.deg.mean);
ActtunDir_loc = vecdir==round(actualnVec.datadegmean);

%if actual tuning not produced by permutations, 
%pick the closest value
if sum(ActtunDir_loc)==0
    [substitutenvecdir,ActtunDir_loc] = min(abs(vecdir - round(actualnVec.datadegmean)));
end

sig = actualnVec.tunStrngth > tunStrngCIbyDir(ActtunDir_loc,2);

%backup
actualnVecStrng = actualnVec.tunStrngth;
