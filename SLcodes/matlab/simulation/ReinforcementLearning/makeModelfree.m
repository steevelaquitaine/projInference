

% Model-free action valuation at each node of
% a tree design(breadth: two-levels)

clear; close all
numTrials=10000;

%Q-learning
QmfS0a1=0.5*ones(numTrials,1);
QmfS0a2=0.5*ones(numTrials,1);
QmfS11a1=0.5*ones(numTrials,1);
QmfS11a2=0.5*ones(numTrials,1);
QmfS12a1=0.5*ones(numTrials,1);
QmfS12a2=0.5*ones(numTrials,1);

PimfS0a1=0.5*ones(numTrials,1);
PimfS11a1=0.5*ones(numTrials,1);
PimfS12a1=0.5*ones(numTrials,1);

S0a1=nan(numTrials,1);
S0a2=nan(numTrials,1);
S11a1=nan(numTrials,1);
S11a2=nan(numTrials,1);
S12a1=nan(numTrials,1);
S12a2=nan(numTrials,1);

R=nan(numTrials,1);

alpha=0.01;
beta=10;

for t=1:numTrials
    %node I
    %action preference
    PimfS0a1(t)=1/(1+exp(-beta*(QmfS0a1(t)-QmfS0a2(t))));
    %action selection
    S0a1(t)=rand<PimfS0a1(t);
    S0a2(t)=S0a1(t)==0;

    %node II
    %action preference
    PimfS11a1(t)=1/(1+exp(-beta*(QmfS11a1(t)-QmfS11a2(t))));
    PimfS12a1(t)=1/(1+exp(-beta*(QmfS12a1(t)-QmfS12a2(t))));
    %action selection
    S11a1(t)=S0a1(t)*(rand<PimfS11a1(t));
    S11a2(t)=S0a1(t)*(S11a1(t)==0);
    S12a1(t)=S0a2(t)*(rand<PimfS12a1(t));
    S12a2(t)=S0a2(t)*(S12a1(t)==0);
    
    %reward collected
    R(t)=...
        S11a1(t)*(rand>0.2)+...
        S11a2(t)*(rand>0.4)+...
        S12a1(t)*(rand>0.7)+...
        S12a2(t)*(rand>0.9);
                                                                                                 
    %action value node I
    QmfS0a1(t+1)=QmfS0a1(t)+ ...
        S0a1(t)*alpha*(S11a1(t)*(QmfS11a1(t)-QmfS0a1(t))+S11a2(t)*(QmfS11a2(t)-QmfS0a1(t)));
    QmfS0a2(t+1)=QmfS0a2(t)+ ...
        S0a2(t)*alpha*(S12a1(t)*(QmfS12a1(t)-QmfS0a2(t))+S12a2(t)*(QmfS12a2(t)-QmfS0a2(t)));
 
    %action value node II
    QmfS11a1(t+1)=QmfS11a1(t)+S11a1(t)*alpha*(R(t)-QmfS11a1(t));
    QmfS11a2(t+1)=QmfS11a2(t)+S11a2(t)*alpha*(R(t)-QmfS11a2(t));
    QmfS12a1(t+1)=QmfS12a1(t)+S12a1(t)*alpha*(R(t)-QmfS12a1(t));
    QmfS12a2(t+1)=QmfS12a2(t)+S12a2(t)*alpha*(R(t)-QmfS12a2(t));
end
fig1=figure;
subp1=subplot(1,2,1);
hold all
hlines(1)=plot(QmfS0a1,':','color',[122/255 16.6/255 228/255],'DisplayName','QmfS0a1');
hlines(2)=plot(QmfS0a2,':','color',[27/255 79/255 53/255],'DisplayName','QmfS0a2');
ylim([0 1])
ylim([0,1])
xlim([0,numTrials])

subp2=subplot(1,2,2);
hold all
hlines(3)=plot(QmfS11a1,'b','DisplayName','QmfS11a1');
hlines(4)=plot(QmfS11a2,'r','DisplayName','QmfS11a2');
hlines(5)=plot(QmfS12a1,'k','DisplayName','QmfS12a1');
hlines(6)=plot(QmfS12a2,'g','DisplayName','QmfS12a2');
legend('show')
hlines(7)=plot(repmat(0.8,numTrials,1),'b');
hlines(8)=plot(repmat(0.6,numTrials,1),'r');
hlines(9)=plot(repmat(0.3,numTrials,1),'k');
hlines(10)=plot(repmat(0.1,numTrials,1),'g');
ylim([0,1])
xlim([0,numTrials])

%behavior
fig2=figure;
subp3=subplot(1,2,1);
plot(PimfS0a1,':','color',[122/255 16.6/255 228/255],'DisplayName','PimfS0a1')
ylim([0 1])
subp4=subplot(1,2,2);
hold all
plot(PimfS11a1,'b:','DisplayName','PimfS11a1')
plot(PimfS12a1,'k:','DisplayName','PimfS11a1')
ylim([0 1])
legend('show')

% model-based
PS0S11=0.5*ones(numTrials,1);
PS0S12=0.5*ones(numTrials,1);
QmbS0a1=0.5*ones(numTrials,1);
QmbS0a2=0.5*ones(numTrials,1);
PimbS0a1=0.5*ones(numTrials,1);

for t=1:numTrials
    %node I
    %action preference
    PimbS0a1(t)=1/(1+exp(-beta*(QmbS0a1(t)-QmbS0a2(t))));
    %action selection
    S0a1(t)=rand<PimbS0a1(t);
    S0a2(t)=S0a1(t)==0;
  
    % state transitions (sample mean)
    PS0S11(t+1)=PS0S11(t)+ S0a1(t)*alpha*(S0a1(t)-PS0S11(t));
    PS0S12(t+1)=PS0S12(t)+ S0a2(t)*alpha*(S0a2(t)-PS0S12(t));

    % action valuation (sample mean)
    QmbS0a1(t+1)=PS0S11(t)*max(QmfS11a1(t),QmfS11a2(t));
    QmbS0a2(t+1)=PS0S12(t)*max(QmfS12a1(t),QmfS12a2(t));
end

figure(fig1)
subplot(subp1)
hold on
plot(QmbS0a1,'color',[122/255 16.6/255 228/255],'DisplayName','QmbS0a1')
plot(QmbS0a2,'color',[27/255 79/255 53/255],'DisplayName','QmbS0a2')
legend('show')
plot(repmat(0.7,numTrials,1))
plot(repmat(0.2,numTrials,1))
ylim([0,1])
xlim([0,numTrials])


% behavior
figure(fig2)
subplot(subp3)
hold on
plot(PimbS0a1,'color',[122/255 16.6/255 228/255],'DisplayName','PimbS0a1')
ylim([0 1])
legend('show')


% because valuation is based on a simulation process: actual mental choice
% (hence max(QmfS11a1(t),QmfS11a2(t))) of all possible trial combinations 
% (in stochastic not all; we use estimates and not real R; we approximate 
% over groups of trials depending on the same action).

% One could formulate "QmbS0a1(t+1)=PS0S11(t)*max(QmfS11a1(t),QmfS11a2(t))"
% as follows:
% I choose S0a1 
% I then estimate PS0S11 
% I choose S11a1
% I estimate QmfS11a1 (end of a mentally simulated trial)
% I then choose S0a1 (again)
% I then estimate PS0S11 (again)
% I choose S11a2
% I estimate QmfS11a2 (end of the second mentally simulated trial)
% I assign maximum proximal action value to distal action.
