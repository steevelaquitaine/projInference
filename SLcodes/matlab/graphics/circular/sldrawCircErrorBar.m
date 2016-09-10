%sldrawCircErrorBar.m
%
%
%
% author: steeve laquitaine
%purpose: circular plot mean vector length and error bar by directions
%
%  usage: 
%
%       
%       nullVecDirection = [0 45 90 135 180 45 180 0];
%       coor = SLpolar2cartesian(nullVecDirection',[1 2 3 4 5 6 7 8],'polar');
%       [vecdir,vecdirPdf,veclenMeanByDir,veclenSemByDir,veclenCIBydir] = slgetVectorStats(coor);
% 
%       sldrawCircErrorBar(vecdir,veclenMeanByDir,veclenCIBydir,maxplot)
% 
%


function sldrawCircErrorBar(vecdir,veclenMeanByDir,veclenCIBydir,maxplot)

%radians
vecdirrad = SLde2r(vecdir,1);

%plot
set(gcf,'color','w')

%set max plot
polar(0,maxplot);

%remove polar outer line
%sldeletePolarOuterLine

for i = 1 : length(vecdir)
    
    %vector direction
    hold on; h = polar(vecdirrad(i),veclenMeanByDir(i),'.r');
    set(h,'MarkerSize',14,'color','k')
    
    %plot confidence interval
    hold on; h = polar([vecdirrad(i) vecdirrad(i)],[veclenCIBydir(i,1) veclenCIBydir(i,2)],'k');
    set(h,'color',[.6 .6 .6])
end




