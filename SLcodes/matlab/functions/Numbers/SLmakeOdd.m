

%SLmakeOdd.m
%
%  author: steeve laquitaine
%    date: 2015/04/22
% purpose: create vector of odd number within xmin xmax range

function x = SLmakeOdd(Xmin,Xmax)

x = [0:1:Xmax]*2 - 1;
x = x(x>=Xmin & x<=Xmax);