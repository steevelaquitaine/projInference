
%SLcircle.m
%
%       $Id: SLcircle.m 750 2012-09-08 07:12:46Z steeve $
%     usage: h=SLcircle(x,y,r,colorc)
%        by: steeve laquitaine
%      date: 140526 last modification 140528
%   purpose: calculate the norm of a vector: pythagoras formula


function h=SLcircle(x,y,r,colorc)

%call for help
if ieNotDefined('x')
    help SLcircle
    return
end

%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
h=plot(x+xp,y+yp,'color',colorc,...
    'linesmoothing','on',...
    'linewidth',0.5);
end