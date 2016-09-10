
%slregOutVarEffectOnInstance.m
%
%
% author: steeve laquitaine
%   date: 151018
%purpose: regress-out the weight of a variable from instances for classification
%
%  usage: 
%
%      %get var1 instances
%      [svol,var1] = getStimvol(v,'myRandomCoh','taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
%      i = getInstances(v,c.myROIinfo,svol,'n=inf',c.startLag,c.blockLen);
%      i1 = i1{1}.classify.instances;
%      stimVol1 = i{1}.classify.instanceVol;
%
%      %get var2 instances
%      [svol,var2] = getStimvol(v,'myRandomDir','taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
%      i = getInstances(v,c.myROIinfo,svol,'n=inf',c.startLag,c.blockLen);
%      i2 = i{1}.classify.instances;
%      stimVol2 = i{1}.classify.instanceVol;
%
%      %regressout var2
%      [instancesregout,var1,stimVol1] = slregOutVarEffectOnInstance(i1,var1,stimVol1,i2,var2,stimVol2);
%
%
%to do: need to carefully check the code again

function [instancesregout,var1,stimVol1] = slregOutVarEffectOnInstance(inst1,var1,stimVol1,inst2,var2,stimVol2)

%make predictor matrix
[instances,predictorM,varsnum,varschar,stimvols] = ...
        slmakeGLMfitPredictorMatrix(inst1,var1,stimVol1,inst2,var2,stimVol2);
    
%regressed-out weights of variable 2 from instances
nClassVar1 = length(var1);
for vox = 1 : size(instances,2) 
    
    %get the weight of each variable
    %glmfit (simple linear regression)
    [weightsByVox{vox},devByvox{vox},statsByVox{vox}] = ...
        glmfit(predictorM,instances(:,vox),'normal','link','identity','constant','on');
    
    %recalculate instance by summing all variables weights (var1, cte and noise) except the ones
    %regressed-out. The first weight is the constant. The
    %subsequent weights are weights from variable 1 classes then from variable 2 classes
    %inst = w1 * x1 + cte + error;
    instancesregoutAll(:,vox) = predictorM(:,1:nClassVar1)*weightsByVox{vox}(2:2+nClassVar1-1) ...
        + weightsByVox{vox}(1) + statsByVox{vox}.resid;
end

%rearrange instances by var1
for iClass = 1 : nClassVar1
    [~,pos] = intersect(stimvols,stimVol1{iClass});
    instancesregout{iClass} = instancesregoutAll(pos,:); 
end