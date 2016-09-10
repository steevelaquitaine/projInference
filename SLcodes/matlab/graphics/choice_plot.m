% steeve laquitaine last version 01/ 10/2008 

close all; %clear all;
clc;


%motor success trials 

s=data(:,2);    % motor success and failure
i1=find(s==1);  % motor success num of trials
f=data(i1,:);   % all success data


% choice trials 

t1=f(:,6);                  % data for the target 1
t2=f(:,8);                  % data for the target 2
i2=find(t1~=t2);             % index of Choice trials (target1 'not=' target2)
choicetrials=f(i2,:);     % choice trials

data=choicetrials;

A=data(:,4);   % action
R=data(:,3);   % reward


% computation of observed learning  

length_filter=15;               %odd number of trial to choose
%filter=ones(1,length_filter);   %for a square filter, non gaussian

sigma=2.5;                                                %variance of the gaussian ?
filter_tmp=[-(length_filter-1)/2:(length_filter-1)/2];  %for a local gaussian running bin of sigma trials
filter=exp(-filter_tmp.^2/(2*sigma^2));                 %for a local gaussian running bin of sigma trials
filter=filter/sum(filter);      %1/1 then 1/2 then 1/3 ...divide by the number of trials

A_for_filter=A(end)*ones(length(A)+length_filter-1);     % Add 'bin' trials to A to obtain A_for_filter.These trials have the same value than the last A value
A_for_filter(1:(length_filter-1)/2)=A(1);                % first A= 'bin' first A_for_filter
A_for_filter((length_filter-1)/2+1:(length_filter-1)/2+length(A))=A; 

p_est=zeros(length(A),1);

for i1=1:length(A);
     p_est(i1)=sum(filter.*A_for_filter([i1:i1+length_filter-1]));
end


plot(p_est,'r'); % Observed learning

