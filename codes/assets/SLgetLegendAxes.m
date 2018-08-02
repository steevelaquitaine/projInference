

%SLgetLegendAxes.m
%
%       $Id: SLgetLegendAxes.m $
%        by: steeve laquitaine
%      date: 140804
%   purpose: get figures' legend axes
%     usage: h = SLgetLegendAxes

function h = SLgetLegendAxes
h=findobj(gcf,'Type','axes','Tag','legend');
