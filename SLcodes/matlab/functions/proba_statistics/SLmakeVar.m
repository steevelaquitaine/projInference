
%SLmakeVar.m
%
%       $Id: SLmakeVar.m $
%        by: steeve laquitaine
%      date: 141106
%   purpose: calculate variance of a probability density function
%            of x
%
%     usage:
%
%           pdf = SLGaussianPdfs(-100:1:100,0,1);
%           varX = SLmakeVar(-100:1:100,pdf)
%
%     note:

function varX = SLmakeVar(x,pdf)

%columns
x = SLmakeColumn(x);
pdf = SLmakeColumn(pdf);

%mean
MeanX = SLmakeMean(x,pdf);

%var
varX = pdf'*((x - MeanX).^2);

