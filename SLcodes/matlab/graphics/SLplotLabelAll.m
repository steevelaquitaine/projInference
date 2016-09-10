
%SLplotLabelAll.m
%
%     author: steeve laquitaine
%       date: 150218
%    purpose: labels all y and x plot axes
%
%      usage:
%
%           SLplotLabelAll('myY','myX')  

function SLplotLabelAll(Myylabel,Myxlabel)

%get axes
ax = SLgetAxes;
ax = sort(ax);

%num rows and col axes
[ncols,nrows] = SLplotNumRowAxes;

%order axes in rows and cols
ax_ordered = reshape(ax,ncols,nrows);
ax_ordered = ax_ordered';

%xlabel
for i = 1 : size(ax_ordered,2)
    xlabel(ax_ordered(end,i),Myxlabel)
end

%ylabel
for i = 1 : size(ax_ordered,1)
    ylabel(ax_ordered(i,1),Myylabel)
end