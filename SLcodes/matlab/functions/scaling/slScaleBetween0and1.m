
%slScaleBetween0and1.m
%
%
% author: steeve laquitaine
%   date: 160204
%purpose: scale x between 0 and 1
%  usage:
%
%       scaledx = slScaleBetween0and1(1:20)
%
%
function scaledx = slScaleBetween0and1(x)


scaledx = (x - min(x(:)))/(max(x(:)) - min(x(:)));

