

%slfmriGetVoxTuningVMfitfmins.m
%
%
%author : steeve laquitaine
%purpose: fit von Mises function to voxel mean responses to a stimulus
%         feature (e.g., motion directions)
%
%usage :
%
%       [u,k,A,minsqe,st] = slfmriGetVoxTuningVMfitfmins(y,var0Val)
%
%       note: scaling all parameter to the unit doesn't produce better fit
%             changing TolFun and TolX either
%             (test with one voxel)
%             initialized k0 at 1
%             initialized k0 at 4

function [u,k,A,minsqe] = slfmriGetVoxTuningVMfitfmins(y,var0Val)

%start amplitude search at mean data amplitude
%initialize parameters with visually checked parameters
A0 = mean(y);
k0 = 4;
u0 = var0Val(y==max(y));

%fit
options = optimset('MaxFunEvals',10000,'MaxIter',100000,'display','iter','TolX',1e-4,'TolFun',1e-4);
[p(1,:),sqe] = fminsearch(@(p) ...
    slgetVMfitsqe(var0Val,y,p),[u0;k0;A0],options);

%====================================[================
%Case there's no preference at all: compare vm fit
%with a simple flat line fit k=0 and A=mean(data))
%and keep the best fit
%====================================================
%note that 180 deg is later set to 'NaN'
p(2,:) = [180 0 mean(y)];
sqe(2) = slgetVMfitsqe(var0Val,y,p(2,:));

%get best parameter set that minimizes
%sum of squared erro between model predictions
%and data
[~,bestpos] = min(sqe);
u = p(bestpos,1);
k = p(bestpos,2);
A = p(bestpos,3);
minsqe = sqe(bestpos);

%plot
% x2rad = SLde2r([1:1:360],0); u2rad = SLde2r(u,0);  
% hold on; plot(A*exp(k.*cos(x2rad-u2rad)-k),'r')

