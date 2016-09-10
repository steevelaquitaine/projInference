

% slInitBayesianModel.m
%
%
%    author: steeve laquitaine
%      date: 160603
%   purpose: setup Bayesian models for model fitting of motion direction
%            estimation data
%
%
%====================
%Basic Bayesian model 
%====================
%fit trial-data with maximum likelihood ('MaxLikelihoodFit')
%
%
%     SLfitBayesianModel({'sub01'},...
%     [80 40 20 1.74 4.77 10.74 34.25 NaN 0.001 15 NaN],...
%     'experiment','vonMisesPrior',...
%     'filename','datafit','MAPReadout','MaxLikelihoodFit','fminsearch');
%
%
%==================================
%Basic Bayesian model with no prior 
%==================================
%Fit basic Bayesian model assuming no prior : with prior strengths k 
%fixed at 0.
%
%
%      SLfitBayesianModel({'sub01'},...
%     [80 40 20 1.74 4.77 10.74 34.25 NaN 0.001 15 NaN],...
%     'experiment','vonMisesPrior',...
%     'filename','datafit','MAPReadout','MaxLikelihoodFit','fminsearch');
%
%
%=======================
%Bayesian sampling model 
%=======================
%Fit Bayesian model that samples posterior
%
%
%     SLfitBayesianModel({'sub01'},...
%     [80 40 20 1.74 4.77 10.74 34.25 NaN 0.001 15 NaN],...
%     'experiment','vonMisesPrior',...
%     'filename','datafit','SamplingReadout','MaxLikelihoodFit','fminsearch');
%