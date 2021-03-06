
%SLmakePredictionsBayesianModel.m
%
% author: steeve laquitaine
%   date: 141122 updated 150110
%purpose: make Bayesian model motion estimates predictions.
%
%usage:
%
%       SLmakePredictionsBayesianModel(displ,coh,pstd,fitP,priorShape,priorModes,TheModel,TrialOrMean,output,varargin)
%

%predictions (trial & average)
function [meanPred, stdPred, cond, PestimateGivenModelUniq,MAP,output] = ...
    SLmakePredictionsBayesianModel(myData,displ,...
    coh,...
    pstd,...
    fitP,...
    priorShape,...
    priorModes,...
    TheModel,...
    TrialOrMean,...
    output,...
    varargin)

[negLogl,~,Logl_pertrial,PestimateGivenModel,MAP,AIC,PdataGivenModel] = SLgetLoglBayesianModel(myData,...
    displ,coh,pstd,fitP,priorShape,...
    priorModes,TheModel,varargin{:});

%output
output.negLogl = negLogl;
output.Logl_pertrial = Logl_pertrial;
output.nTrials = length(displ);
output.AIC = AIC;
output.PdataGivenModel = PdataGivenModel;

%case von Mises prior
%--------------------
%conditions
if strcmp(priorShape,'vonMisesPrior')
    [cond,idxCondUniqtrial,~] = SLuniqpair([pstd coh displ]);
end

PestimateGivenModelUniq = PestimateGivenModel(:,idxCondUniqtrial);

%sort everything to match data sorting (ascending order)
[cond,PosSorted] = sortrows(cond,[-2 -1]);
PestimateGivenModelUniq = PestimateGivenModelUniq(:,PosSorted);

%predictions about estimate mean and std (circular mean and std) for
%each task condition (columns)
%-------------------------------------------------------------------
%Use circular statistics
numCond = size(PestimateGivenModelUniq,2);
meanPred = nan(numCond,1);
stdPred = nan(numCond ,1);

for i = 1 : numCond
    data = SLcircWeightedMeanStd(MAP, PestimateGivenModelUniq(:,i));
    meanPred(i)= data.deg.mean;
    stdPred(i) = data.deg.std;
end

%case trial-predictions
%----------------------
if strcmp(TrialOrMean,'Trial')

    %case von Mises Prior
    %--------------------
    if strcmp(priorShape,'vonMisesPrior')
        output.TrialPred = nan (size(pstd));
        for i = 1 : numCond
            trialsThiC = pstd==cond(i,1) & coh==cond(i,2) & displ==cond(i,3);
            output.TrialPred(trialsThiC) = meanPred(i);
        end
    end
end
