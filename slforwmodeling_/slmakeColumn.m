
%SLmakeColumn.m
%
% author: steeve laquitaine
%   date: 1407725
%purpose: make column vector
%
%  usage: x = SLmakeColumn(x)

function x=SLmakeColumn(x)

if size(x,1)<size(x,2)
    x=x';
end
