
% SLerrorarea.m
%
% author: steeve laquitaine 
%   date: 12/03/2010 start: 17:30 ; end :18:17 - Debug1 16.03.2010; end 11:44
%purpose: plot error area around mean
%
%  usage:
%
%       h = SLerrorarea([1:10],linspace(3,3,10),1:10,'input');
%
%       when 'input', arg1 = mean or median, arg2 = error, arg3 = x.   

function h = SLerrorarea(arg1,arg2,arg3,arg4,arg5,varargin)

%case default
if sum(strcmp(varargin,'default'))
    theMeanOrMedian = nanmean(arg1,2);
    theMeanOrMedian(isnan(theMeanOrMedian)==1) = [];
    
    %sem
    for i = 1 : size(arg1,1)
        error(i,1) = sem(arg1(i,:));
    end
    
    error(inan) = [];
    Eup = theMeanOrMedian + error;
    Edn = theMeanOrMedian - error;
end

%case mean, and error (e.g., std or sem..) are input
if sum(strcmp(varargin,'input'))
    theMeanOrMedian = SLmakeColumn(arg1);
    error = SLmakeColumn(arg2);
    Eup = theMeanOrMedian + error;
    Edn = theMeanOrMedian - error;
    x = arg3;
end

%case mean, and error up and error down are input
if sum(strcmp(varargin,'inputErrors'))
    theMeanOrMedian = SLmakeColumn(arg1);
    error = SLmakeColumn(arg2);
    Eup = SLmakeColumn(arg4);
    Edn = SLmakeColumn(arg5);
    x = arg3;
end
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





