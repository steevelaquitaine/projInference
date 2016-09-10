
%SLcopyPlotData.m
%
%
%author: steeve laquitaine
%  date: 150605
%purpose: copy x and y data of a plot
%
%  usage :
%
%           select a plot (arrow)
%           [xdata,xdata] = SLcopyPlotData

function [xdata,ydata] = SLcopyPlotData

xdata=get(gco, 'xdata'); 
ydata=get(gco,'ydata');
