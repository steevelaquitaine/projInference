

%sem.m
%
%author: steeve laquitaine 
%  date: 12042009


function y=sem(x,idx)



%UNTITLED1 Summary of this function goes here
%  Detailed explanation goes here

%  use a matrix
%  output :Standard error of the mean std/sqrt(count) for each line of the matrix


%sem over rows
if idx==1
    for i=1: size(x,1)
        x2=x(i,isnan(x(i,:))==0);
        sp_sz=size(x2,2);
        y(i)=nanstd(x2)/sqrt(sp_sz-1);%sem
    end
else
    %sem over col
    for i=1: size(x,2)
        x2=x(isnan(x(:,i))==0,i);
        sp_sz=size(x2,1);
        y(i)=nanstd(x2)/sqrt(sp_sz-1);%sem
    end
end
