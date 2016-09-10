%slfmriGetSingleVoxNullvmfitTuningDist.m
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
%          1000 permutations of 100 valued-y takes 40 secs
%
%
%usage:
%
%          nVec = slfmriGetSingleVoxNullvmfitTuningDist(rand(10,1),...
%               [0 90 180 270 0 90 180 270 0 90]',1000)
%
%y and x must be column vectors

function [nVec,actualnVec,minsqe,A] = slfmriGetSingleVoxNullvmfitTuningDist(y,x,nperm)

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
    [ymean_by_xu{i},xu] = slCalculateMeanSortedBy1Var(yperm,x);         
    
    %von mises fit    
    [datadegmean(i),tunStrngth(i),A(i),~] = slfmriGetVoxTuningVMfitfmins(ymean_by_xu{i},xu);
end


%% to visualize fits
% i = 10
% hold on; plot(xu,ymean_by_xu{i},'o')
% u = datadegmean(i); k = tunStrngth(i); a = A(i);
% x2rad = SLde2r([1:1:360],0); u2rad = SLde2r(u,0);  
% hold on; plot(a*exp(k.*cos(x2rad-u2rad)-k),'r')
%%

%backup
%null tunings
nVec.deg.mean = round(datadegmean);
nVec.deg.tunStrngth = tunStrngth;

%voxel tuning
[ymean_by_xu,xu] = slCalculateMeanSortedBy1Var(y,x);
[actualnVec.datadegmean,actualnVec.tunStrngth,A,minsqe] = slfmriGetVoxTuningVMfitfmins(ymean_by_xu,xu);
actualnVec.A = A;
