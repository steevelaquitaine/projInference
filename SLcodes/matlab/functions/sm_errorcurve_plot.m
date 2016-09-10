function p=sm_errorcurve_plot(xaxis,x,xx,xcolor)
% steeve last version 12/ 10/2008 

% use a matrix
% give upper sem and lower sem for error curve (same principle as errorbar)

% xx row or col


% get ou nan
x=x(1:max(max(size(x,1)-sum(isnan(x),1))),:);
%mean
MEAN=nanmean(x',xx)';
%SEM
for i=1:size(x,1)
    SEM(i,1)=nanstd(x(i,:))/sqrt((length(x(i,:))-(sum(isnan(x(i,:)))))-1);
end
%up
up= MEAN+SEM;
down= MEAN-SEM;
p=[up MEAN down]
p=sgolayfilt(p,1,11);
plot(xaxis,p,xcolor)
