%slfmriGetSingleVoxNullgridSearchfitTuningDist.m
%
% author:  steeve laquitaine
%purpose:  calculate single vox null neural vector distribution, by
%          permutation of responses y association to stimulus circular 
%          feature x (e.g., motion direction, orientation etc...) and voxel
%          actual neural vector.
%
%          The number of unique permutations of y vis a vis x is too large
%          for long vector y (e.g., y of 100 values, number of permutations 
%          of )
%
%usage:
%
%          nVec = slfmriGetSingleVoxNullgridSearchfitTuningDist(rand(10,1),...
%               [0 90 180 270 0 90 180 270 0 90]',1000)
%
%y and x must be column vectors

function [nVec,actualnVec,minsqe] = slfmriGetSingleVoxNullgridSearchfitTuningDist(y,x,nperm)

num = length(y);

%preallocate for speed
datadegmean = nan(nperm,1);

%get voxel tuning and strength for each 
%response permutation
parfor i = 1 : nperm
    
    %permute vox responses
    j = randperm(num);
    yperm = y(j);   
    
    %get mean y responses for each x  
    %Averaging removes the effect of the non-uniform 
    %distribution of x on the voxel feature vector (0.001)        
    [ymean_by_xu,xu] = slCalculateMeanSortedBy1Var(yperm,x);         
    
    %von mises fit    
    [datadegmean(i),tunStrngth(i),~,~] = slfmriGetVoxTuningVMfitGridSearch(ymean_by_xu,xu);
end

%backup
%null tunings
nVec.deg.mean = round(datadegmean);
nVec.deg.tunStrngth = tunStrngth;

%voxel tuning
[ymean_by_xu,xu] = slCalculateMeanSortedBy1Var(y,x);
[actualnVec.datadegmean,actualnVec.tunStrngth,actualnVec.A,minsqe] = slfmriGetVoxTuningVMfitGridSearch(ymean_by_xu,xu);

