 
%SLcircShiftXandY4Plot.m.m
%
% author: steeve laquitaine
%   date: 1407725
%purpose: recenter by circular shift x and y vectors.
%
%  usage: 
%         [xcentered,ycentered,x4plot] = SLcircShiftXandY4Plot(1:10,1:10,-5)
%         plot(x4plot,ycentered)
%         set(gca,'xtick',x4plot,'xticklabel',xcentered)


function [xshifted,yshifted,x4plot] = SLcircShiftXandY4Plot(x,y,shift)

%make sure x and y are column vectors
x=SLmakeColumn(x);
y=SLmakeColumn(y);

%x for plot
x4plot=SLmakeColumn(1:1:numel(x));

%x and y centered
xandycentered=circshift([x y],[shift,0]);
xshifted=xandycentered(:,1);
yshifted=xandycentered(:,2);
