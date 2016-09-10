

%SLpolar2cartesian.m
%
%   Date: 131007  last modified: 150722
% Author: Steeve Laquitaine
%Purpose: Convert angle (in deg or rads) to cartesian coordinates (x,y)
%
%
%Inputs:
%
%     theta :  angle
%         r :  unit circle radius
%
% Usage:
%
%      coord = SLpolar2cartesian([1:1:360]',1)
%
%
%varargin:
%
%           'polar'
%           'radian'


function [coord,theta2] = SLpolar2cartesian(theta,r,varargin)

%check args
if nargin <2    
    fprintf('\n (SLpolar2cartesian) Please indicate the radius    \n')
    keyboard    
end

%convert polar to radian
if sum(strcmp(varargin,'polar')) || isempty(varargin)
    
    theta = SLmakeColumn(theta);    %column vector
    theta2.deg = theta;             %angle    
    theta2.rad = single(theta2.deg)*pi/180; %from degree to radian
    
    %verbose
    if sum(strcmp(varargin,'verbose'))

        fprintf('(SLpolar2cartesian) Assuming input angles are in deg.  \n')
        fprintf('(SLpolar2cartesian) Theta converted to rads then to carts. \n')
    
    end

elseif sum(strcmp(varargin,'radian'))
    
    theta2.deg = SLde2r(theta,0);   %degrees
    theta2.rad = theta;             %radians
    
    %verbose
    if sum(strcmp(varargin,'verbose'))
        fprintf('(SLpolar2cartesian) Theta is in radian and is converted to cartesian.  \n')
    end
end

%Convert radians to cartesian
theta2.rad = SLmakeColumn(theta2.rad);
r = SLmakeColumn(r);
x = r.*cos(theta2.rad);
y = r.*sin(theta2.rad);

%round x and y to 10-4 to prevent weird angle due to imprecision
x = round(x*10^4)/10^4;
y = round(y*10^4)/10^4;
coord = [x y];



