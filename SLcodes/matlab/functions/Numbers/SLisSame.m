
%SLisSame.m
%
% author: Steeve Laquitaine
%   date: 140425
%purpose: find if two vectors are same

function o=SLisSame(x1,x2)

x1=SLmakeColumn(x1);
x2=SLmakeColumn(x2);

%check size and content
if numel(x1) == numel(x2)
    o = sum(x1==x2)==numel(x1);
else
    o = 0;
end



