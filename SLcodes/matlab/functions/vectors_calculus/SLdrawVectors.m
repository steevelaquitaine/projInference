% drawVectors.m
%
% author: steeve laquitaine
%   date: 131128 last modif: 140602
%
%  usage: coor=polar2cartesian([180 90],[1 1])
% "theta" is an angle in degree
% "r" is the radius of the unit circle (vector magnitude)
%
% reference: http://www.mathsisfun.com/polar-cartesian-coordinates.html

function drawVectors(theta,Vlength)
polar(0,1)

%draw vectors
coor=polar2cartesian(theta,Vlength);
arrow([0,0],[coor(1) coor(2)]);