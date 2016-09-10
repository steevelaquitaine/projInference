

%slgetAllpolarOjects.m
%
%
% author: steeve laquitaine
%purpose: get axis object handles for texts, lines, axes and patch
%
%usage: 
%
%         [texts, lines, axess, patches] = slgetAllpolarOjects

function [texts, lines, axess, patches] = slgetAllpolarOjects

%get axis object handles for texts, lines, axes and patch
texts = findall(gca,'type','text');
lines = findall(gca,'type','line');
axess = findall(gca,'type','axes');
patches = findall(gca,'type','patch');
