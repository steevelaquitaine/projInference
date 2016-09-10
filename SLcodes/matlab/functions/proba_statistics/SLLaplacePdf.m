

% SLLaplacePdf.m
%
%     author: Steeve Laquitaine
%       date: 141223
%    purpose: create a Laplace probability distribution
%
%      usage:
%
%           p = SLLaplacePdf(1:1:400,225,1,'norm')
%
%
%Description: 
%
%       p = (1./2.*b).*exp(- abs(x - mu)./b)
%
% inputs: 'norm' scales to probabilities, otherwise [];


function p = SLLaplacePdf(x,mu,b,type)

%Laplace
p = (1./2.*b).*exp(- abs(x - mu)./b);

%scale to probabilities
if strcmp(type,'norm')
    p = p./sum(p);
end

p = SLmakeColumn(p);