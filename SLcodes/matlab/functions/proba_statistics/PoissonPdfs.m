
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

poispdf = poisspdf(x,lamda);

%scale to pdfs.
if strcmp(type,'norm')==1
    Z_ = sum(poispdf);
    Z = Z_(ones(numel(x),1),:);
    poispdf  = poispdf ./Z;
else
end