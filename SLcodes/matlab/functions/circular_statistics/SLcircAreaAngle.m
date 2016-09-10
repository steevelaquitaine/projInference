
%SLcircAreaAngle.m
%
%
% author: steeve laquitaine
%   date: 150723
%purpose: highlight angle (area) between two vectors
%
%  usage:
%
%           SLcircAreaAngle(0,90)
%           SLcircAreaAngle(0,90,'EdgeColor','none','FaceAlpha',0.1)


function SLcircAreaAngle(theta1,theta2,mycolor,varargin)

set(gcf,'color','w')

%cartesian
th1coor = SLpolar2cartesian(theta1,1);
th2coor = SLpolar2cartesian(theta2,1);

%note: patch builds from y(1) to y(end) back to y(0)
x = [0 th1coor(1) th2coor(1) 0];
y = [0 th1coor(2) th2coor(2) 0];
patch(x,y,mycolor,varargin{:})
axis off
xlim([-1 1])
ylim([-1 1])