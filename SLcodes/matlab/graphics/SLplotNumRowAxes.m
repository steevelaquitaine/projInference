
%SLplotNumRowAxes.m
%
%     author: steeve laquitaine
%       date: 150218
%    purpose: get number of row and col axes in a plot
%
%usage
%
%   SLplotNumRowAxes

function [nCols, nRows] = SLplotNumRowAxes

%get axes
ax = SLgetAxes;

%remove legends handles
%----------------------
lgH = SLgetLegendAxes;
ax = setdiff(ax,lgH);

%case nb axes > 1
if length(ax) > 1
    pos = cell2mat(get(ax,'position'));
else
    %case 1 axis
    pos = get(ax,'position');
end
nCols = numel(unique(pos(:,1)));
nRows = numel(unique(pos(:,2))); 

