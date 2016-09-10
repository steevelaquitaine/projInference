
%slGetAngleYAsclosestToX.m
%
%
% author: Steeve Laquitaine
%purpose: Transform angles y (deg) in the linear space so that the are 
%         expressed as the closest angles to x (in signed deg)
%         e.g., y = 355 and x = 10 are very distant on the linear space. We
%         can rewrite y as -5 so that it is closer to x in the linear
%         space.
%
%  usage: 
%
%       newY = slGetAngleYAsclosestToX(355,10)        
%
%input:
%       Y : scalar or vector to be converted to a new angle value 
%       X : scalar or vector reference for converting Y


function newY = slGetAngleYAsclosestToX(x,y)

%transform y in the linear space so that it is expressed as the closest 
%angle to x
shortestDist = SLvectors2signedAngle(y,x,'polar');
newY = x + shortestDist;


