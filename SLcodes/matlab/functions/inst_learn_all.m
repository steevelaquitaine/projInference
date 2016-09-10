close all; clear all;clc

load A1;
load R1;


A=A1;
R=R1;

t0s=find(A==0&R==1);
t0f=find(A==0&R==0);
t1s=find(A==1&R==1);
t1f=find(A==1&R==0);

dt=1;

exlude_trials=[length(A):-1:length(A)-dt+1];

t0s=setdiff(t0s,exlude_trials);
t0f=setdiff(t0f,exlude_trials);
t1s=setdiff(t1s,exlude_trials);
t1f=setdiff(t1f,exlude_trials);

%number of suboptimal choice rewarded/nonrewarded followed by
%optimal/suboptimal choice
n0s0=sum(A(t0s+dt)==0); % when A0 is rewarded how much is it chosen on 'dt' next trial 
n0s1=sum(A(t0s+dt)==1); % when A0 is rewarded how much does it switch to A1 on the 'dt' next trial
n0f0=sum(A(t0f+dt)==0); % when A0 is not rewarded how much is it chosen on 'dt' next trial
n0f1=sum(A(t0f+dt)==1); % when A0 is not rewarded how much does it switch to A1 on the 'dt' next trial

[n0s0 n0f0;n0s1 n0f1]

[n0s0/(n0s0+n0s1) n0f0/(n0f0+n0f1)] % if reward driven choice,  n0s0/(n0s0+n0s1) > n0f0/(n0f0+n0f1) and [n0s0 and n0f1] > [n0f0 and n0f0]
n0s0_f0=[n0s0/(n0s0+n0s1) n0f0/(n0f0+n0f1)];

%number of optimal choice rewarded/nonrewarded followed by
%optimal/suboptimal choice
n1s0=sum(A(t1s+dt)==0);
n1s1=sum(A(t1s+dt)==1);
n1f0=sum(A(t1f+dt)==0);
n1f1=sum(A(t1f+dt)==1);

[n1s0 n1f0;n1s1 n1f1]

[n1s0/(n1s0+n1s1) n1f0/(n1f0+n1f1)] % should be very low
n1s0_f0=[n1s0/(n1s0+n1s1) n1f0/(n1f0+n1f1)];






%% Create axes
axes1 = axes('OuterPosition',[0 0.4985 0.4823 0.5015],'XTick',[1 2 3 4],'XtickLabel',{'n0s0','n0f0','n0s1','n0f1'}); 
title(axes1,'After A0');
hold(axes1,'all');
%% Create bar
bar1 = bar(1,n0s0,'Parent',axes1,'FaceColor',[1 0 0]);
bar2 = bar(2,n0f0,'Parent',axes1,'FaceColor',[0 0 0]);
bar3 = bar(3,n0s1,'Parent',axes1,'FaceColor',[1 0 0]);
bar4 = bar(4,n0f1,'Parent',axes1,'FaceColor',[0 0 0]);
%% Create axes
axes2 = axes('OuterPosition',[0 0 0.4823 0.4985],'XTick',[1 2 3 4],'XtickLabel',{'n1s0','n1f0','n1s1','n1f1'}); 
title(axes2,'After A1');
hold(axes2,'all');
bar5 = bar(1,n1s0,'Parent',axes2,'FaceColor',[1 0 0]);
bar6 = bar(2,n1f0,'Parent',axes2,'FaceColor',[0 0 0]);
bar7 = bar(3,n1s1,'Parent',axes2,'FaceColor',[1 0 0]);
bar8 = bar(4,n1f1,'Parent',axes2,'FaceColor',[0 0 0]);
%% Create axes
axes3 = axes('OuterPosition',[0.4823 0.4985 0.5177 0.5015],'XTick',[1 2],'XtickLabel',{'n0s0/(n0s0+n0s1)','n0f0/(n0f0+n0f1)'}); 
title(axes3,'Stay A0 after S vs F');
hold(axes3,'all');
bar9 = bar(1,n0s0/(n0s0+n0s1),'Parent',axes3,'FaceColor',[1 0 0]);
bar10 = bar(2,n0f0/(n0f0+n0f1),'Parent',axes3,'FaceColor',[0 0 0]);
%% Create axes
axes4 = axes('OuterPosition',[0.4823 0 0.5177 0.4985],'XTick',[1 2],'XtickLabel',{'n1s0/(n1s0+n1s1)','n1f0/(n1f0+n1f1)'}); 
title(axes4,'Switch to A0 after S vs F');
hold(axes4,'all');
bar11 = bar(1,n1s0/(n1s0+n1s1),'Parent',axes4,'FaceColor',[1 0 0]);
bar12 = bar(2,n1f0/(n1f0+n1f1),'Parent',axes4,'FaceColor',[0 0 0]);
legend('rewarded','omitted')
set(gcf,'Color',[0.9 0.9 1])
