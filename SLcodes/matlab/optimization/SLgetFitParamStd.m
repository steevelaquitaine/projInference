
%SLgetFitParamStd.m
%
% author: Steeve Laquitaine
%   date: 140829
%purpose: get the std of a model's parameters
%
%  usage:
%
%       fitParamStd = SLgetFitParamStd(Hessian,data,pred,fitP)
%
%Hessian can be obtained with optimization algorithm such as "lsqnonlin" or
%"fmincon". If you have the jacobian instead you can use it too to calculate
%the covariance matrix of the fit parameters: d.covar = inv(jacobian'*jacobian);
%
%note: Std with complexe values means parameters have negative variance

function fitParamStd = SLgetFitParamStd(Hessian,data,pred,fitP)

%Covariance matrix of the parameters
%d.covar = inv(jacobian'*jacobian);
d.covar = inv(Hessian); 

%noise variance
residual = (data - pred)';
residual = residual(~isnan(residual));
numP = sum(~isnan(fitP));
numData = sum(~isnan(data)); 
noiseVariance = (residual*residual')/(numData - numP);

%std of model parameters
fitParamStd = diag(sqrt(noiseVariance*d.covar))';