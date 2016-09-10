close all; clear all;clc

load action004010

num_trials=size(action,2);
last_trials=[];
num_av=50;

for i1=1:num_trials
    i2=min(find(isnan(action(:,i1))==1));
    if i2>num_av
        last_trials=[last_trials; (action(i2-num_av:i2-1,i1))'];
    end
end

coldel=find(size(action,1)-sum(isnan(action))<= num_av);
colkeep=find(size(action,1)-sum(isnan(action))>= num_av)';
Mean_lasts=mean(last_trials,2);

table=[colkeep Mean_lasts];

%plot(mean(last_trials,2),'o')

num_bins=10;
epsi=1e-5;
dx=(1+2*epsi)/num_bins;
x_range=-epsi+dx*[0:num_bins];

h=histc(mean(last_trials,2),x_range);

bar(x_range(1:end-1)+dx/2,h(1:end-1))

p=mean(mean(last_trials,2));



p_expect=-ones(1,num_av+1);

for i1=0:num_av
    p_expect(i1+1)=nchoosek(num_av,i1)*p^i1*(1-p)^(num_av-i1);
end

% for i1=1:num_bins;
%     i2=find

hold on
plot([0:num_av]/num_av,p_expect*num_av*num_trials/num_bins,'r')

Xlabel('Average performance on the last 50 trials');
Ylabel('Number of session')