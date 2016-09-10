
%sldrawVecPdf.m
%
% author: steeve laquitaine
%purpose: plot vector probability distribution in a circular graph
%
%  usage: 
%
%       
%       nullVecDirection = [0 45 90 135 180 45 180 0];
%       coor = SLpolar2cartesian(nullVecDirection',[1 2 3 4 5 6 7 8],'polar');
%       [vecdir,vecdirPdf] = slgetVectorStats(coor);
% 
%       sldrawVecPdf(vecdir,vecdirPdf,maxplot)
% 
function sldrawVecPdf(vecdir,vecdirPdf)

%radians
vecdirrad = SLde2r(vecdir,1);

%plot
set(gcf,'color','w')
mp = max(vecdirPdf);
polar(0,mp);

%remove polar outer line
sldeletePolarOuterLine

%vector direction
for i = 1 : length(vecdir)
    hold on; h(i) = polar([0 vecdirrad(i)],[0 vecdirPdf(i)],'-k');
end
set(h,'linewidth',3,'color',[.6 .6 .6])