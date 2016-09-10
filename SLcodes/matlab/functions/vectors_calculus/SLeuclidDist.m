
%SLeuclidDist.m
%
% author: Steeve Laquitaine
%   date: 141229
%purpose: Euclidean distance between two data points
% 
%  usage:  
%
%       e = SLeuclidDist(rand(10,1),rand(10,1))
%
%options (varargin)

function e = SLeuclidDist(v1,v2)

%col vectors
v1 = SLmakeColumn(v1);
v2 = SLmakeColumn(v2);

%Euclidean distance
%based on Pythagorean theorem
e = sqrt(sum((v1 - v2).^2));











