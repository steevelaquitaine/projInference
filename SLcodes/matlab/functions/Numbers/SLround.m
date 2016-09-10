
%SLround.m
%
% author: steeve laquitaine
%   date: 150128
%purpose: round to floating unit
%
%  usage:
%
%           data = SLround(5.87,0.1)

function data = SLround(data,Floating)

data = round(data*(1./Floating))./(1./Floating);



