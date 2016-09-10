
%Fit the data with a linear model
function [slope,intercept] = linefit(x,y,xnew,linecolor)

%available data
i=~isnan(y);
y=y(i);

%fit
xfit = x;
if size(x,1)~=size(y,1)
    x = x';
end
[slope,intercept] = polyfit(x(i),y,1);
yfit = polyval(slope,xnew);
plot(xnew,yfit,'-',...
    'linewidth',1,...
    'color',linecolor,...
    'LineSmoothing','on');