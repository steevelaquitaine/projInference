
%slsubstRowAndcolVec.m
%
%
% author : steeve laquitaine
%   date : 160905
%purpose : quickly substract row and col vectors x and y
%          in the following way
%
%          x = 1:10; y = x';
%          for i = 1 : length(y)
%             xminy(i,:) = x - y(i);
%          end
%
%usage :
%
%       xminy = slsubstRowAndcolVec(1:10,[1:10]')


function xminy = slsubstRowAndcolVec(x,y)

x = SLmakeRow(x);
y = SLmakeColumn(y);

xminy = bsxfun(@(x,y) x-y,s,y);