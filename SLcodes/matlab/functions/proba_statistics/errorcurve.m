function p=errorcurve(x)
% steeve last version 12/ 10/2008 

% use a matrix
% give upper sem and lower sem for error curve (same principle as errorbar)

clc;

%mean
for i=1:length(x(:,1))
    MEAN(i,1)=nanmean(x(i,:))
end

%SEM
for i=1:length(x(:,1))
    sem(i,1)=nanstd(x(i,:))/sqrt((length(x(1,:))-(sum(isnan(x(i,:)))))-1);
end
    
%up
up= MEAN+sem;
down= MEAN-sem;

p=[up down]

