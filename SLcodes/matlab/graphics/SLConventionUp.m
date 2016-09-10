
%SLConventionUp.m
%
% author: Steeve Laquitaine
%   date: 140804
%purpose: clean up all figures' axes with the lab conventions. 
%         e.g., drawpublishAxis, etc...
%
%   usage: 
%
%       fig = figure(1);
%       SLConventionUp

function SLConventionUp

%select figure
figure(gcf)

%get figure's axes
h = SLgetAxes;

%get rid of legend axes
hlg = SLgetLegendAxes;
h = setdiff(h,hlg);

%convention up each axis
numAxes = length(h);
for i = 1 : numAxes
    axes(h(i))
    SLdrawPublishAxis
end

SLremoveDeadSpace