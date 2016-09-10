
%SLmakeColumn.m
%
% author: steeve laquitaine
%   date: 1407725
%purpose: recenter by circular shift x and y vectors.
%
%  usage: [xcentered,ycentered,x4plot] = SLmakeColumn(x)

function SLmakeColumn(x)

if size(x,1)<size(x,2)
    x=x';
end
