

%slLinearizeCircDataForScatterAndLinFit.m
%
% author: steeve laquitaine
%   date:
%purpose: linearize y circular data as a function of x circular data
%  usage:
%
%       [x,y] = slLinearizeCircDataForScatterAndLinFit(225,270)
%
%Description:
%
%     Calculate "d" the shortest distance between "y" and "x"
%     y = d + x


function [x,y] = slLinearizeCircDataForScatterAndLinFit(x,y)

%Transform "y" as "x" + the shortest distance vector to "x"
y = x + slvectors2signedAngle(y,x,'polar');



