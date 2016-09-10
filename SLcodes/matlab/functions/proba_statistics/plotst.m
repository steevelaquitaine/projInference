
% steeve last version 12/ 10/2008 

%close all; %clear all;
clc;


% computation of observed learning  
A=data;

sigma=5;
range_filter=length(A);
filter_tmp=(-range_filter:range_filter)';
filter=exp(-filter_tmp.^2/(2*sigma^2));

p=-ones(1,length(A));
for i1=1:length(A);
    p(i1)=nansum(A.*filter(range_filter-i1+2:2*range_filter-i1+1))./nansum(filter(range_filter-i1+2:2*range_filter-i1+1));
end

p=p';

%plot(p,'-+')    

%hold on

%plot(A,'r.')

clear('A','sigma','range_filter','filter_tmp','filter','i1');
