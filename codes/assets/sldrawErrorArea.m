
% sldrawErrorArea.m
%
%
% author: steeve laquitaine
%   date: 14/12/2016 
%purpose: plot error area around mean
%
%  usage:
%
%       h = sldrawErrorArea([1:10],linspace(3,3,10),1:10,'input');
%
%   inputs :
%
%
%           Case we have a vector of means and a vector of error:
%               y: a vector of means
%               e: a vector of error values
%               x: a vector of x-axis values
%           varargin : 'color', [1 0 0 ]


function h = sldrawErrorArea(y,e,x,varargin)

theMeanOrMedian = SLmakeColumn(y);
error = SLmakeColumn(e);
Eup = theMeanOrMedian + error;
Edn = theMeanOrMedian - error;
x = SLmakeColumn(x);

%graphics
%color mean
if sum(strcmp(varargin,'color'))
    pos = find(strcmp(varargin,'color'))+1;
    meancolorOrig = varargin{pos};
else
    %default
    meancolorOrig = [.7 .7 .7];
end

%darken mean and whitened area
meancolor = SLdarken(meancolorOrig,0.5);
fillcolor = SLwhiten(meancolorOrig,0.7);

%----
%plot
%----
%area
hold all;
fill([x; x(end:-1:1)],[Eup; Edn(end:-1:1)],...
    fillcolor,'edgecolor',fillcolor)

%mean
h = plot(x,theMeanOrMedian,'color',meancolor,...
    'linesmoothing','on','linewidth',1);





