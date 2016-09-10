
%SLmakeRow.m
%
% author: steeve laquitaine
%   date: 140904
%purpose: make row vector
%
%  usage: x = SLmakeRow(x)

function x=SLmakeRow(x)

if size(x,1)>size(x,2)
    x=x';
end
