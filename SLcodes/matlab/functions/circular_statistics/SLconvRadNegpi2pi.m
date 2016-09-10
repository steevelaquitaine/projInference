
%SLconvRadNegpi2pi.m
%
% author: steeve laquitaine
%   date: 19/06/2015
%purpose: convert radians to values ranging between -pi and pi.
%
%  usage: 
%
%       [rad,input] = SLconvRadNegpi2pi(-pi:0.1:2*pi)

function [rad,input] = SLconvRadNegpi2pi(rad)

%initialize
rad2 = rad;

%warning
if any(rad2(:) > 2*pi) || any(rad2(:) < -pi)
    fprintf('(SLconvRadNegpi2pi) The min/max input value of a radian is -pi to 2*pi. Your values are out of that range \n')
    keyboard
end

%convert
rad2(rad > pi) = rad(rad > pi) - 2*pi;
rad2(rad < - pi) =  rad(rad < - pi) + 2*pi;

input = rad;
rad = rad2;
