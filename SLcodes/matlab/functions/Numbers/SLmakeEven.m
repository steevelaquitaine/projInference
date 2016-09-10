

%SLmakeEven.m
%
%  author: steeve laquitaine
%    date: 2015/04/22
% purpose: create vector of even number within xmin xmax range

function x = SLmakeEven(Xmin,Xmax)

x = [0:1:Xmax]*2;
x = x(x>=Xmin & x<=Xmax);
x(x==0)=[];