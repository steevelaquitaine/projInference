

%sllinefit.m
%
% author: Steeve Laquitaine
%purpose: Fit the data with a linear model
%
%varargin:
%   
%   'dispoff': no display

function [slope,intercept] = sllinefit(x,y,xnew,linecolor,varargin)

linewidth = 1;
if any(strcmp(varargin,'linewidth'))==1
    linewidth = varargin{find(strcmp(varargin,'linewidth'))+1};
end

%available data
i=~isnan(y);
y=y(i);

%fit
xfit = x;
if size(x,1)~=size(y,1)
    x = x';
end
%get slope and intercept of best linear fit
[slope,intercept] = polyfit(x(i),y,1);
%predicted y
yfit = polyval(slope,xnew);
%case plot
if ~any(strcmp(varargin,'dispOff'))
    plot(xnew,yfit,'-','linewidth',linewidth,'color',linecolor,'LineSmoothing','on');
end