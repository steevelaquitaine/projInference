



function SLgraphCopyFig2Fig

h = gcf;
hh = gca;
c = copyobj(hh,get(hh,'Parent'));

i = gcf;

set(c,'Parent',i);

%figure 1 - axis 1 - object 1
xdata = get(gco, 'xdata'); ydata =get(gco,'ydata');

%figure 2  - axis 1 - object 1
hold all; plot(xdata,ydata,'g') 