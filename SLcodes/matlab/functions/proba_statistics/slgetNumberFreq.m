

%slgetNumberFreq.m
%
% author: steeve laquitaine
%   date: 160509
%purpose: get number frequencies
%
%  usage:
%
%     f = slgetNumberFreq([1 1 1 2 2 3 4])

function f = slgetNumberFreq(x)
u = unique(x);
u = u(~isnan(u));
f = [u,histc(x(:),u)];