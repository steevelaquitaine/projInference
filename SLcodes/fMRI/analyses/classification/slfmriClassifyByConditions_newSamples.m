
%slfmriClassifyByConditions_newSamples.m
%
%
% author: steeve laquitaine
%purpose: classify motion directions from the brain activity for data
%         sorted by a condition (e.g., "switching") to determine whether the
%         sensory evidence is represented or not when the subjects switch
%         to his prior.
%         The # of instances used for classification is equalized
%         between the conditions to allow a fair comparison of the
%         classification accuracy between the conditions.
%
%  usage :
%
%           slfmriInitClassifAnalysisTaskDotDirfMRI05
%           params = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances',...
%                     'leaveOneOut','fisher','balancByRemovI=1',...
%                     'numInstances',8,'loadInstancesByCond'}
%           myConds = {{'myRandomDir_x_mySwitch=1_x_myRandomCoh=0.06'},....
%                       {'myRandomDir_x_mySwitch=2_x_myRandomCoh=0.06'}};           
%           [ci,m,s,acc,accSTE,cbyc,sbyc,o]= slfmriClassifyByConditions_newSamples(myConds,params,100,o)
% 
%no random seed is used


function [ci,m,s,acc,accSTE,cbyc,sbyc,o] = slfmriClassifyByConditions_newSamples(myConds,params,nIter,o)

t0 = tic;

%store accuracies and STE for each condition
for i = 1 : nIter    
    [acc{i},accSTE{i},cbyc{i},sbyc{i},o] = slfmriClassifyByConditions(o,myConds,params);
    accAll(i,:) = cell2mat(acc{i});
end

%get mean accuracy and confidence interval over samples of the instances
parfor i = 1 : size(accAll,2)
    [err_margin(i),ci{i},m(i),s(i)] = slMakeCI(accAll(:,i),0.95);  
end

duration = toc(t0);
fprintf('%s %.02f \n','Duration : ',duration)


