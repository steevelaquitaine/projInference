
%slTextFindExp.m
%
%author: steeve
%
% usage:
%   o = slTextFindExp({'Thank you','Maybe not'},'Maybe')


function o = slTextFindExp(text,exp)

o = regexp(text,exp);