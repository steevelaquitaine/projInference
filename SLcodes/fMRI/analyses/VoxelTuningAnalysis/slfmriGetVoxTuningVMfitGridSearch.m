

%slfmriGetVoxTuningVMfitGridSearch.m
%
%
%author : steeve laquitaine
%purpose: fit von Mises function to voxel mean responses to a stimulus
%         feature (e.g., motion directions)
%         The resolution of the tuning selectivity is 5 deg to speed 
%         up the code (4x faster than 1 deg)
%
%usage :
%
%       [u,k,A,minsqe,st] = slfmriGetVoxTuningVMfitGridSearch(y,var0Val)


function [u,k,A,minsqe] = slfmriGetVoxTuningVMfitGridSearch(y,var0Val)

%start amplitude search at mean data amplitude
%k at 4 and selectivity at 180
%'TolFun',1e-1,'TolX',1e-1 are the best speed/fit quality tradeoff
%parameters. They produces good qualitative results within
%10 min/500 voxel (parallel computing)
options = optimset('TolFun',1e-1,'TolX',1e-1,'display','on');
A = nanmean(y);

%steps of 5 degrees to speed up the code
us = 1:5:360;
parfor i = 1:1:length(us)
    [p(i,:),sqe(i)] = fminsearch(@(p) ...
        slgetVMfitsqe(var0Val,y,p),[us(i);4;A],options);
end
%====================================================
%Case there's no preference at all: compare vm fit
%with a simple flat line fit k=0 and A=mean(data))
%and keep the best fit
%====================================================
%note that 180 deg is later set to 'NaN'
p(length(us)+1,:) = [180 0 mean(y)];
sqe(length(us)+1) = slgetVMfitsqe(var0Val,y,p(length(us)+1,:));

%get best parameter set that minimizes
%sum of squared error between model predictions
%and data
[minsqe,bestpos] = min(sqe);
u = p(bestpos,1);
k = p(bestpos,2);
A = p(bestpos,3);
minsqe = sqe(bestpos);