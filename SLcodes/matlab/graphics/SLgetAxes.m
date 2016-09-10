

%SLgetAxes.m
%
%       $Id: SLgetAxes.m $
%        by: steeve laquitaine
%      date: 140709
%   purpose: get figures axes
%     usage: 
%           myAxes = SLgetAxes

function myAxes = SLgetAxes
myAxes = findall(gcf,'type','axes');
