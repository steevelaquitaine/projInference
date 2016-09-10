

%SLfindBetween.m
%
%author: steeve laquitaine
%purpose: find position of values of x and values of x between 
%         "BinStart" and "BinEnd"
% update: 160503
%
%  usage:
%
%       pos = SLfindBetween(10:-1:1,0,2,'includeBinEnd')
%       pos = SLfindBetween(1:1:10,2,4,'includeBinEnd')

function [pos,vals] = SLfindBetween(x,BinStart,BinEnd,type)

if strcmp(type,'includeBoth')
    pos = intersect(find(x >= BinStart),find(x <= BinEnd));
    vals = x(pos);
end
if strcmp(type,'includeBinStart')
    pos = intersect(find(x >= BinStart),find(x < BinEnd));
    vals = x(pos);
end
if strcmp(type,'includeBinEnd')
    pos = intersect(find(x > BinStart),find(x <= BinEnd));
    vals = x(pos);
end
if strcmp(type,'includeNone')
    pos = intersect(find(x > BinStart),find(x < BinEnd));
    vals = x(pos);
end

