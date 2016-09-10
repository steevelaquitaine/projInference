%slfmriClassifyBootInstByConditions.m
%
%
% author: steeve laquitaine
%purpose: classify motion directions from the brain activity for data
%         sorted by a condition (e.g., "switching") to determine whether the 
%         sensory evidence is represented or not when the subjects switch 
%         to his prior.
%         The # of instances used for classification is equalized
%         between the conditions to allow a fair comparison of the
%         classification accuracy between the conditions by sampling n trials
%         over the instances many times and averaging over their decoding
%         accuracy.
%
%
%setup
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%
%           params = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances',...
%               'leaveOneOut','fisher','balancByRemovI=1',...
%               'numInstances',8,'loadInstancesByCond'};
%
%           myConds = {{'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'},....
%               {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'}};
%
%           [ci1,m1,s1,ci2,m2,s2,accAll,accSTEAll,cbyc,sbyc] = slfmriClassifyBootInstByConditions(100,params,myConds,o)


function [ci1,m1,s1,ci2,m2,s2,accAll,accSTEAll,cbyc,sbyc] = slfmriClassifyBootInstByConditions(numBoot,params,myConds,o)

%bootstrap and classify
parfor i = 1 : numBoot
    [acc,accSTE,cbyc{i},sbyc{i}] = slfmriClassifyByConditions(o,myConds,params);
    accAll(i,:) = cell2mat(acc);
    accSTEAll(i,:) = cell2mat(accSTE);
end

%Get accuracy stats
[err_margin1,ci1,m1,s1] = slMakeCI(accAll(:,1),0.95);
[err_margin2,ci2,m2,s2] = slMakeCI(accAll(:,2),0.95);
