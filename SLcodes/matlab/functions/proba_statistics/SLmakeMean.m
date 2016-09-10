
%SLmakeMean.m
%
%       $Id: SLmakeMean.m $
%        by: steeve laquitaine
%      date: 141106
%   purpose: calculate standard deviation of a probability density function
%            of x
%
%     usage:
%           x = -100:1:100;
%           p = SLGaussianPdfs(x,0,1);
%           meanX = SLmakeMean(x,p)

function meanX = SLmakeMean(x,pdf)

%column
x = SLmakeColumn(x);
pdf = SLmakeColumn(pdf);

%mean(projection of x into pdf)
meanX = pdf'*x;