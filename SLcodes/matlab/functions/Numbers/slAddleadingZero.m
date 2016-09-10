

%slAddleadingZero.m
%
%author: steeve laquitaine
%  date: 160420
%purpose : add a leading zero to an number
%
%e.g., x = slAddleadingZero(5)
%
%>> 05

function x = slAddleadingZero(x)

x = sprintf('%02d',x);