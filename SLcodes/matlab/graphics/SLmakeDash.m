
%SLmakeDash.m
%
% author: steeve laquitaine
%   date: 140904
%purpose: make clearer dashed lines
%usage:
%
%       yDashed = SLmakeDash(1:360,1:360)
%
%reference
%http://www.mathworks.com/matlabcentral/answers/57632-long-dashes-in-a-dashed-line-matlab-plot

function yDashed = SLmakeDash(x,y)

%make sure column vectors
x = SLmakeRow(x);
y = SLmakeRow(y);

%start point
a = 1;  
%gap width
b = 2;  
%gap frequency
c = 10; 

rangeX = max(x) - min(x);
rangeY = max(y) - min(y);
Q = [diff(x) / rangeX; diff(y) / rangeY];
L = [0, cumsum(sqrt(sum(Q .* Q)))];
L = (length(L) / L(end)) * L;
index   = rem(round(L) + a, c) <= b;
yDashed = y;
yDashed(index) = NaN;
