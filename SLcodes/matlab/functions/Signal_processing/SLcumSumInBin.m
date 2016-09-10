
% author: Steeve Laquitaine
%   date: 140807
%purpose: Calculate sum within bins.
%
%  usage:
%       p = [1 1 1 2 2 2 3 3 3]';
%       bins = [1 1 1 2 2 2 3 3 3]';
%       cumSumInBin = SLcumSumInBin(y,binsIdx)

function cumSumInBin = SLcumSumInBin(y,binsIdx)

cumSumInBin = accumarray(binsIdx,y,[],@(x) sum(x));



