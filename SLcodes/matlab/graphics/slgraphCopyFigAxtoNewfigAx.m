
%slgraphCopyFigAxtoNewfigAx.m
%
%
% author: steeve laquitaine
%purpose: copy an axis from one figure to another
%
%usage: 
%
%       figure(1);
%       oldAxis = plot(rand(5));
%       figure(2);
%       newAxis = subplot(1,2,1);
%       slgraphCopyFigAxtoNewfigAx(oldAxis,figure(2),newAxis)


function slgraphCopyFigAxtoNewfigAx(oldAxis,newFig,newAxis)

%copy old axis to new figure
c = copyobj(oldAxis,newFig);

%set new axis position
%in new figure
pos = get(newAxis,'position');
delete(newAxis)
set(c,'position',pos)