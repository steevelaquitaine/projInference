

%SLsetFigPublishSize.m
%
%
% author: steeve laquitaine
%   date: 141219
%purpose: set figure size in centimeters for publication print.
%         8.7 cm wide, free height.
%
%reference: 
%
%   http://www.nature.com/srep/authors/submit.html#general-figure

function SLsetFigPublishSize

%get current position
cpos = get(gcf,'position');

%size of figure's "drawing" area on screen
set(gcf, 'Units','centimeters', 'Position',[0 0 8.7 8.7])

%size on printed paper
set(gcf, 'PaperPositionMode','auto')


%set back to original position on screen
%---------------------------------------
%new height and width in pixels
set(gcf, 'Units','pixels')
newpos = get(gcf,'position');
NewHeiht_Width = newpos(3:4);

%set back to original position
set(gcf,'position',[cpos(1:2) NewHeiht_Width])
