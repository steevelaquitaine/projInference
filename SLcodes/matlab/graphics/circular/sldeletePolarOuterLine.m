

%sldeletePolarOuterLine.m
%
%
% author: steeve laquitaine
%purpose: remove the polar thick outerline in a figure
%         get axis object handles for texts, lines, axes and patch
%
%usage: 
%
%       sldeletePolarOuterLine

function sldeletePolarOuterLine

%get axis object handles for texts, lines, axes and patch
a = findall(gca,'type','text');
b = findall(gca,'type','line');
c = findall(gca,'type','axes');
d = findall(gca,'type','patch');

%remove polar outer line
set(d,'visible','off')
% set(b(8),'visible','off')
set(gcf,'color','w')