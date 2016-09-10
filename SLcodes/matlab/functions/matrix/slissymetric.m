
%slissymetric.m
%
%
% author: steeve laquitaine
%purpose: check if a matrix is symmetric
%
%usage: 
%
%       o = slissymetric(eye(10,10))
%       o = slissymetric(rand(10,10))


function o = slissymetric(SIGMA)

o = isequal(SIGMA,SIGMA.');