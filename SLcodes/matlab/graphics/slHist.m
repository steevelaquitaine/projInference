
%slHist.m
%
% author: steeve laquitaine
%   date: 150913
%purpose: customized version of hist
%
%usage:
%
%       [n,x] = slHist(rand(100,1),0:0.2:1)

function [n,x] = slHist(y,x)

hist(y,x);
set(get(gca,'child'),'FaceColor',[.6 .6 .6],'EdgeColor','none');
set(gca,'xtick',x,'xticklabel',x)
xlim([0 1])
box off