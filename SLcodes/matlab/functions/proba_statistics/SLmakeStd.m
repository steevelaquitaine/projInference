
%SLmakeStd.m
%
%       $Id: SLmakeStd.m $
%        by: steeve laquitaine
%      date: 141106
%   purpose: calculate standard deviation of a probability density function
%            of x
%
%     usage:
%
%           pdf = SLGaussianPdfs(-100:1:100,100,1);
%           stdX = SLmakeStd(-100:1:100,pdf)
%
%     note:

function stdX = SLmakeStd(x,pdf)

%columns
x = SLmakeColumn(x);
pdf = SLmakeColumn(pdf);

%mean
varX = SLmakeVar(x,pdf);

%std
stdX = sqrt(varX);