


function makeGabor1D
%A gabor multiply a sinusoid & a gaussian

%space
x=-30:1:30;

%gaussian
m=0;
s=1;
Yg=normpdf(x,m,s);
figure('color','w')
hold all; plot(x,Yg,':','color',[.7 .7 .7],'linewidth',2)

%sinusoid
Ys=sin(x);
plot(x,Ys,'color',[.8 .8 .8],'linewidth',2)

%gabor
Ygabor=Yg.*Ys;
plot(x,Ygabor,'color','k','linewidth',6)

%title
title('1Dgabor')

drawPublishAxis