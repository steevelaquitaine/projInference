
%SLsimCompetitionPredictionsWithExpPrior.m
%
% author: steeve laquitaine
%   date: 141122 updated 150110
%purpose: Simulate Bayesian estimate predictions with actual experimental
%         priors and motion directions used in a motion experiment.
%
%usage:
%
%       SLsimCompetitionPredictionsWithExpPrior([72.5 7.95 1.24 3.7 3.701 inf inf NaN 0.001 16 NaN],...
%             'bimodalPrior')
%
%predictions:
%   - 12% coh (kllh ~ 7.95) best maximize qualitatively the three peaks in estimates distributions
%   - prior modes strength ~ 3.7 best maximizes the three peaks and allows motion
%   directions far from the modes (conflict)

function SLsimCompetitionPredictionsWithExpPrior(simParameters,varargin)

%experimental conditions
priorShape = 'bimodalPrior';
TheModel = 'withoutCardinal';

%simulate experiment with two bimodal priors 
%-------------------------------------------
%bimodal prior 1
modes1 = [145 305];
o = SLinitRunExpBimodalPrior('sim',225,modes1,simParameters(4),[.06 .12 .24],[107 75 33],5:10:355);
d = o.parameter.dir.series';
coh = o.parameter.dir.coh';
pstd = repmat(o.parameter.dir.strength,numel(d),1);
priorModes = repmat(modes1,numel(d),1);

%bimodal prior 2
modes2 = [165 285];
o = SLinitRunExpBimodalPrior('sim',225,modes2,simParameters(5),[.06 .12 .24],[107 75 33],5:10:355);
d = [d ; o.parameter.dir.series'];
pstd = [pstd; repmat(o.parameter.dir.strength,numel(o.parameter.dir.series),1)];
coh  = [coh ; o.parameter.dir.coh'];
priorModes = [priorModes; repmat(modes2,numel(o.parameter.dir.series),1)];
output.uniqCond = SLuniqpair([pstd coh d]);
numcond = size(output.uniqCond,1);


%Predictions
%status
fprintf('%s \n','Simulating model predictions...')

%models' estimate mean, std and distribution predictions
output.fitP = simParameters;
[meanPred,stdPred,cond,PestimateGivenModelUniq,...
    MAP] = SLmakePredictionsCompetitionModel(d,...
    coh,...
    pstd,...
    output.fitP,...
    priorShape,...
    priorModes,[],...
    output,...
    varargin);

%make sure predicted and data distributions are calculated on
%the same space. Works for predicted estimate space = 1:1:360 deg.
%adjust predicted estimate probability distribution from 1:1:360 to
%0:10:360 by summing probabilities within consecutive bins of 10
%deg (law of probabilities).
%commonSpace = 0:10:360;
commonSpace = 0:10:360;
numSpace = numel(commonSpace)-1;
[~,bins] = histc(0:1:360,commonSpace);
bins(end) = [];
bins = bins';
numCond = size(cond,1);
pPred = nan(numSpace,numCond);
for ijk = 1 : numCond
    pPred(:,ijk) = SLcumSumInBin(PestimateGivenModelUniq(:,ijk),bins);
end
commonSpace = commonSpace(2:end);

%predicted mean,std and distributions
output.meanPred    = meanPred;
output.stdPred     = stdPred;
output.meanDisPred = pPred;
output.uniqCond    = cond;

%draw predicted mean, std and distribution
%mean and std
SLdrawModelsPredictionCentered([],[],output.meanPred,output.stdPred,...
    output.uniqCond,priorModes,priorShape)

%distribution
SLdrawModelsPredictionHistwithDataCentered([],commonSpace,...
    output.meanDisPred,...
    d,coh,pstd,priorModes,output.uniqCond,priorShape,varargin)
