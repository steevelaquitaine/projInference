
% SLPoissonPdfs.m
%
%      author: Steeve Laquitaine
%        date: 141105
%       usage: poispdf = SLPoissonPdfs(1:1:360,180,'norm')
%
%      inputs: 'norm' scales to probabilities, otherwise [];
%               when k are same, u must be a value of x for the code to work.
%
% description: lamda is both the mean and variance that are equal for a
%              Poisson distribution.

function poispdf = SLPoissonPdfs(x,lamda,type)

%Poisson distribution is supposed 
%to characterize integer values
if any(SLisinteger(x)==0)
    sprintf('(SLPoissonPdfs) x is not a set of integers')
end

poispdf = poisspdf(x,lamda);
poispdf = SLmakeColumn(poispdf);

%scale to pdfs.
if strcmp(type,'norm')==1
    Z = sum(poispdf);
    poispdf  = poispdf./Z;
else
end