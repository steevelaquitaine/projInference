

%   SLisinteger.m
%
%     author: steeve laquitaine
%       date: 140711
%
%    purpose: state if a number is an integer (1) or not (0)
%
%      usage: output=isinteger(1,2)
%      
%   

function output=SLisinteger(y)

output=y-fix(y)==0;