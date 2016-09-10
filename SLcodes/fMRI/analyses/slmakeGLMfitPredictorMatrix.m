
%slmakeGLMfitPredictorMatrix.m
%
%    author: steeve laquitaine
%      date: 151018
%   purpose: make a predictor matrix for a vector of observations
%decription:
%         
%        instances: Nobservations by Ndimensions (e.g., voxels) responses to predict
%             var1:
%
%  usage:
%      
%      %(fMRI)
%      %get variable 1 instances
%      [svol, var1] = getStimvol(v,'myRandomCoh','taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
%      i = getInstances(v,c.myROIinfo,svol,'n=inf',c.startLag,c.blockLen);
%      i1 = i1{1}.classify.instances;
%      stimVol1 = i{1}.classify.instanceVol;
%
%      %get variable 2 instances
%      [svol, var2] = getStimvol(v,'myRandomDir','taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
%      i = getInstances(v,c.myROIinfo,svol,'n=inf',c.startLag,c.blockLen);
%      i2 = i{1}.classify.instances;
%      stimVol2 = i{1}.classify.instanceVol;
%
%      %prediction matrix
%      [instances,predictornum,predictors,stimvols] = slregOutVarEffectOnInstance(i1,var1,stimVol1,i2,var2,stimVol2)


function [instances,predictorM,varsnum,varschar,stimvols] = slmakeGLMfitPredictorMatrix(inst1,var1,stimVol1,inst2,var2,stimVol2)


%var1
instmp   = [];
vartmp  = [];
svoltmp = [];
for i = 1 : length(inst1)
    %instances
    instmp  = [instmp;inst1{i}];
    %get var (a number)
    var1num = str2num(var1{i}(find(var1{i}=='=')+1:end));
    vartmp  = [vartmp;repmat(var1num,size(inst1{i},1),1)];
    %vols
    svoltmp = [svoltmp stimVol1{i}];  
end

%var2
instmp2  = [];
vartmp2  = [];
svoltmp2 = [];
for i = 1: length(inst2)
    instmp2  = [instmp2 ; inst2{i}];
    var2num  = str2num(var2{i}(find(var2{i}=='=')+1:end));
    vartmp2  = [vartmp2 ; repmat(var2num,size(inst2{i},1),1)];
    svoltmp2 = [svoltmp2 stimVol2{i}];  
end

%match stimvols, var1 and var2
[~,pos1] = sort(svoltmp);
[~,pos2] = sort(svoltmp2);
instmp  = instmp(pos1,:); 
vartmp  = vartmp(pos1,:);
svoltmp = svoltmp(pos1);
instmp2  = instmp2(pos2,:); 
vartmp2  = vartmp2(pos2,:);
svoltmp2 = svoltmp2(pos2);

%create predictor matrix
predictorM = zeros(length(svoltmp),length([var1 var2]));
var1u = unique(vartmp);
for i = 1 : length(var1)
    predictorM(find(vartmp==var1u(i)),i)=1;    
end
var2u = unique(vartmp2);
for i = 1 : length(var2)
    predictorM(find(vartmp2==var2u(i)),length(var1)+i)=1;    
end

%check that all stimvols are matched to two predictors
if any(sum(predictorM,2)~=2)
   fprintf('%s \n',['(slmakeGLMfitPredictorMatrix) WARNING, all ' ...
                    'observation should be matched with exactly two ' ...
                    'predictors !'])    
   keyboard
end

%output for glmfit
varschar  = [var1 var2];
varsnum   = [var1u' var2u'];
instances = instmp;
stimvols  = svoltmp'; 

fprintf('%s \n', '(slmakeGLMfitPredictorMatrix) Predictor matrix...done.')