
%slGetVectorLength.m
%
%
%author : steeve laquitaine
%purpose: calculate the length of a vector from his cartesian coordinates x
%         and y


function vecLen = slGetVectorLength(x,y)

vecLen = sqrt(sum([x y].^2));
