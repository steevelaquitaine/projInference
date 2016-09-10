function p=errorcurve_plot(x,xcolor)
% steeve last version 12/ 10/2008 

% use a matrix
% give upper sem and lower sem for error curve (same principle as errorbar)

% xx =column; always 2
% x=data
% xaxis= x axis
clc;

x1=1;
xx=2;

% get ou nan
%x=x(1:max(max(size(x,x1)-sum(isnan(x),x1))),:);

%mean
MEAN=nanmean(x,xx);

%SEM
for i=1:size(x,x1)
    sem(i,1)=nanstd(x(i,:))/sqrt((length(x(i,:))-(sum(isnan(x(i,:)))))-1)';
end    
sem(find(isnan(MEAN)==1))=NaN;


%up down sem
up= MEAN+sem;
down= MEAN-sem;
p=[up MEAN down];
plot(p,xcolor)
