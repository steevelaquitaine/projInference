function [coord] = polar2cartesian(theta,r)
% theta is an angle in degree
% r is the radius of the unit circle
% Coord are in visual angle
% Record angle in degree
theta2.deg = theta;
% Convert from degree to radian
theta2.rad = theta2.deg*pi/180;
% Calculate visual angles coordinates
x = r*cos(theta2.rad);
y = r*sin(theta2.rad);
coord = [x y];
