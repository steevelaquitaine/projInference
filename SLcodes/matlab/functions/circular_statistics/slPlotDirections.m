
%slPlotDirections.m
%
%
% author: Steeve laquitaine
%   date: 150722
%purpose: Plot directional arrows.
%         theta can be in deg, rads or cartesian
%
%  usage:
%
%       slPlotDirections([1;2;3;4],'polar')
%       slPlotDirections([pi;p./2;],'radian')
%       slPlotDirections([pi;p./2;],'cartesian')
%       slPlotDirections(theta,unit,'polar','Property1',PropVal1,'Property2',PropVal2,...)
%
function slPlotDirections(theta,unit,varargin)

if nargin<2
    fprintf('\n (slPloyDirections) Please indicate the angle unit \n')
    fprintf('(slPloyDirections) You choices are                   \n')
    fprintf('(slPloyDirections) - "polar"                         \n')
    fprintf('(slPloyDirections) - "radian"                        \n')
    keyboard
end

if strcmp(unit,'polar')
    
    set(gcf,'color','w')
    carts = SLpolar2cartesian(theta,1);         %cartesian coordinates
    arrow([0 0],carts,'linesmoothing','on',varargin{:})     %direction
    xlim([-1 1])
    ylim([-1 1])
    axis square
    axis off
    
elseif strcmp(unit,'radian')
    
    set(gcf,'color','w')
    carts = SLpolar2cartesian(theta,1,'radian');  %cartesian coordinates
    arrow([0 0],carts,'linesmoothing','on')       %direction
    xlim([-1 1])
    ylim([-1 1])
    axis square
    axis off
    
elseif strcmp(unit,'cartesian')
    
    set(gcf,'color','w')
    arrow([0 0],theta,'linesmoothing','on') %direction
    xlim([-1 1])
    ylim([-1 1])
    axis square
    axis off
    
end

