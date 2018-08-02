
%SLdarken.m
%
% author: Steeve Laquitaine
%   date: 140904
%purpose: make color darker
%
%  usage:
%
%     aDarkened = SLdarken([0 0 .5],0.2)

function aDarkened = SLdarken(colorC,tune)
aDarkened = (0-colorC).*tune + colorC;