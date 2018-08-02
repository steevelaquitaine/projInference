

% author: Steeve Laquitaine
%   date: 140904
%purpose: make color whiter
%
%  usage:
%
%     aWhitened = SLwhiten([0 0 .5],0.7)

function aWhitened = SLwhiten(colorC,tune)
aWhitened = (1-colorC).*tune + colorC;