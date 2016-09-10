
%author: steeve laquitaine
  %date: 140411
 %usage: SLpositionFigure(fig,0.5,0.5)
     %e.g., fig is figure handle obtain with gcf;
     %MoveHorizontal=0.5 mmove figure on the horizontal axis by 50% of its
     %actual position.
     
function SLpositionFigure(fig,MoveHorizontal,MoveVertical)

%get figure position
position=get(gcf,'Position');

%position somewhere on the screen
%divide the screen in subplot as for the subplot function and position
%figure as in subplot function.

%re-position figure
set(fig,'position',[MoveHorizontal*position(1) MoveVertical*position(2) position(3) position(4)],'color','w')
