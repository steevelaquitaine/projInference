

%slGetAIC.m

% author: steeve laquitaine
%   date: 150802
%purpose: Calculate Akaike information criterion (AIC) for model comparison. 
%		  It is a metric of the likelihood (logl) of data given a model with nP parameters.
%		  The likelihood is penalized by nP such that data are less likely to come from
%		  high nP models. 
%		  The model with the lowest AIC is the most likely model.
%
%  usage:
%
%		AIC = slGetAIC(9,30000)

function AIC = slGetAIC(nP,logl)

AIC = 2 * (nP - logl)