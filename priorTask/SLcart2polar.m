     

%SLcart2polar.m
% 
%      author: steeve laquitaine
%        date: 131128 last modif 140714
%
%     purpose: convert from cartesian coordinates (x,y) to polar coordinates
%               (theta,radius).
%               radius=sqrt(x^2+y^2)
%               rad=atan2(y,x)
%
%       usage: 
%
%               [rad,deg,radius] = SLcart2polar([1 0; 0.5 0.5; 0 1; -0.5 0.5; -1 0; -0.5 -0.5; 0.5 -0.5])
%               x and y can be < 0.


function [rad,deg,radius] = SLcart2polar(v)

%radius
x = v(:,1);
y = v(:,2);
radius = sqrt(x.^2 + y.^2);

%angle (radian & degrees)
%When x or y are negative angle must be adjusted according to his quadrant 
%(in degree) otherwise atan produces wrong values.
rad = atan2(y,x);
deg = SLra2d(rad);
