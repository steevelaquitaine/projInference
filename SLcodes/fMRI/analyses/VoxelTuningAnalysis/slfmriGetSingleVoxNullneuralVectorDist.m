
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
%          nVec = slfmriGetSingleVoxNullneuralVectorDist([10 1 4 7],[0 90 180 270],1000)

function [nVec,actualnVec] = slfmriGetSingleVoxNullneuralVectorDist(y,x,nperm)

%calculate neural vector for shuffled responses
num = length(y);

%preallocate for speed
datadegmean = nan(nperm,1);
datacoordmean = nan(nperm,2);

parfor i = 1 : nperm
    
    %permute vox responses
    j = randperm(num);
    yperm = y(j);   
    
    %average responses by x  
    %we need to average responses to remove the effect
    %of the non-uniform distribution of x on the voxel
    %feature vector (0.001)        
    [ymean_by_xu,xu] = slCalculateMeanSortedBy1Var(yperm,x);     
    [datadegmean(i),datacoordmean(i,:)] = slfmriCalculateVoxFeatureVector(xu,ymean_by_xu/nansum(ymean_by_xu));
end


% tic
% %permute responses
% perIdx = cell2mat(arrayfun(@(dummy) randperm(num),1:nperm,'UniformOutput',false)');
% yperm = y(perIdx);
% 
% %locate each condition
% [xu,~,cond_location] = unique(x);
% 
% %average data for each condition
% %meanbyCond = accumarray(cond_location,yperm,[length(conds) 1],@mean);
% %get conditions row and col indices
% meanbyCond = slCalcMeanByRowCondEachCol(yperm,cond_location);
% p = bsxfun(@rdivide,meanbyCond,nansum(meanbyCond,1));
% 
% %preallocate for speed
% datadegmean = nan(nperm,1);
% datacoordmean = nan(nperm,2);
% parfor i = 1 : nperm
%     [datadegmean(i),datacoordmean(i,:)] = slfmriCalculateVoxFeatureVector(xu,p(:,i));
% end
% toc

nVec.deg.mean = round(datadegmean);
nVec.coord.mean = datacoordmean;

%voxel actual neural vector
[ymean_by_xu,xu] = slCalculateMeanSortedBy1Var(y,x);
[actualnVec.deg.mean,actualnVec.coord.mean] = slfmriCalculateVoxFeatureVector(xu,ymean_by_xu/nansum(ymean_by_xu));

