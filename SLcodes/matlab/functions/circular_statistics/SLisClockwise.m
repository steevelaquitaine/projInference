
%SLisClockwise.m
%
% author: Steeve laquitaine
%   date: 150722
%purpose: get if x is clockwise to y
%
%  usage: 
%
%       clockwise = SLisClockwise(90,45,'polar')
%       clockwise = SLisClockwise(0,pi./2,'radian')
%
%x should be a scalar or column vector
%y should be a scalar

function clockwise = SLisClockwise(x,y,type)

if nargin < 3
    
    fprintf('(SLisClockwise) Please indicate angle type. \n')
    fprintf('(SLisClockwise) Your choices are: \n')
    fprintf('(SLisClockwise) - radian \n')
    fprintf('(SLisClockwise) - polar \n')
    
    keyboard
    
end

if strcmp(type,'radian')

    x = SLra2d(x);
    y = SLra2d(y); %convert to deg
    
elseif strcmp(type,'polar')
   
end

angle = SLvectors2signedAngle(x,y,'polar');  %get angle distance
clockwise = angle<0; %clockwise ?
