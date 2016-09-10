%SLmakeR2.m
%
%
% author: steeve laquitaine
%   date: 140830
%purpose: calculate the R-squared which is the percent of variance in a 
%         dataset that is explained by a model
%
% usage:
%
%   R2 = SLmakeR2(data,SSE)
%
%inputs: 
%      data: the data used to fit the model 
%       SSE: sum of squared error between the data and the model
%            predictions
%
%outputs:
%       R2: the R-squared
%
%
%Description:
%   R-squared = 1 - SSE/SST
%   SSE (alias residual) is the SSE calculated for the best model's prediction.
%   note: R^ (alias coefficient of determination) can be <0 when the model
%   does an awful job predicting the data. In this case SSE exceeds that is
%   the model fits the data even worse than does a horizontal line.
%   Model may not be appropriate or constraints may not be set correctly.
%   see http://www.graphpad.com/support/faqid/711/
%
%Another way is:
%   function R2 = makeR22(data,pred)
%   %R2=1 - SSE/SST
%   %SSE (alias residual) is the SSE calculated for the best model's prediction.
%   %note: R^ (alias coefficient of determination) can be <0 when the model
%   %does an awful job predicting the data. In this case SSE exceeds that is
%   %the model fits the data even worse than does a horizontal line.
%   %Model may not be appropriate or constraints may not be set correctly.
%   %see http://www.graphpad.com/support/faqid/711/    
%

%model's R-square
function R2 = SLmakeR2(data,SSE)

%R=corr(data,pred);
%R2=R^2;
SST = sum((data - nanmean(data)).^2);
R2 = 1 - SSE/SST;

%case perfect fit
if SSE==0
    R2 = 1;
end