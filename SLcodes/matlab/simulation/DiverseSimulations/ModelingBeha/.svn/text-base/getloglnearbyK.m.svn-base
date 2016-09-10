

%% Look at the -sum(log likelihood) of parameters neighboring the true 
%parameters. The data are simulated data which we know the true parameters.
%We expect to see min(-sum(log likelihood)) lying at the true parameters.
tic
%kl1=350;
upo=upo';
kl1=0:50:700;
logL=nan(numel(kl1),1);
for i=1:numel(kl1);
    [logL(i)]=makeSSE5(upo,d,coh,pstd,[kl1(i) 24 7.5 4.3 8.9 14.5 33]);
end
figure('color','w')
hold all
scatter(kl1,logL,100,'ko','markerfacecolor','w','linewidth',1)
[minlogL,I]=min(logL);
plot([350 350],[0 logL(kl1==350)],':','color',[.7 .7 .7],'linewidth',2)
text(350,-10000,'350','color',[.4 .4 .4])
plot([kl1(I) kl1(I)],[0 minlogL],'r:','linewidth',2)
text(kl1(I),-10000,num2str(kl1(I)),'color','r')
title('-Sum(log(likelihood)) nearby kl1(350)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kl2=24
kl2=0:24:700;
logL2=nan(numel(kl2),1);
for i=1:numel(kl2);
    [logL2(i)]=makeSSE5(upo,d,coh,pstd,[350 kl2(i) 7.5 4.3 8.9 14.5 33]);
end
figure('color','w')
hold all
scatter(kl2,logL2,100,'ko','markerfacecolor','w','linewidth',1)
[minlogL,I]=min(logL2);
plot([24 24],[0 logL2(kl2==24)],':','color',[.7 .7 .7],'linewidth',2)
text(24,-10000,'24','color',[.4 .4 .4])
plot([kl2(I) kl2(I)],[0 minlogL],'r:','linewidth',2)
text(kl2(I),-10000,num2str(kl2(I)),'color','r')
title('-Sum(log(likelihood)) nearby kl2(24)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kl3=7.5
kl3=0:7.5:700;
logL3=nan(numel(kl3),1);
for i=1:numel(kl3);
    logL3(i)=makeSSE5(upo,d,coh,pstd,[350 24 kl3(i) 4.3 8.9 14.5 33]);
end
figure('color','w')
hold all
scatter(kl3,logL3,100,'ko','markerfacecolor','w','linewidth',1)
plot([7.5 7.5],[0 logL3(kl3==7.5)],':','color',[.7 .7 .7],'linewidth',2)
text(7.5,-10000,'7.5','color',[.4 .4 .4])
[minlogL,I]=min(logL3);
plot([kl3(I) kl3(I)],[0 minlogL],'r:','linewidth',2)
text(kl3(I),-10000,num2str(kl3(I)),'color','r')
title('-Sum(log(likelihood)) nearby kl3(7.5)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kp1=4.3
hold all
kp1=0:4.3:700;
logL4=nan(numel(kp1),1);
for i=1:numel(kp1);
    [logL4(i)]=makeSSE5(upo,d,coh,pstd,[350 24 7.5 kp1(i) 8.9 14.5 33]);
end
figure('color','w')
hold all
scatter(kp1,logL4,100,'ko','markerfacecolor','w','linewidth',1)
plot([4.3 4.3],[0 logL4(kp1==4.3)],':','color',[.7 .7 .7],'linewidth',2)
text(4.3,-10000,'4.3','color',[.4 .4 .4])
[minlogL,I]=min(logL4);
plot([kp1(I) kp1(I)],[0 minlogL],'r:','linewidth',2)
text(kp1(I),-10000,num2str(kp1(I)),'color','r')
title('-Sum(log(likelihood)) nearby kp1(4.3)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kp2=8.9
kp2=0:8.9:700;
logL5=nan(numel(kp2),1);
for i=1:numel(kp2);
    [logL5(i)]=makeSSE5(upo,d,coh,pstd,[350 24 7.5 4.3 kp2(i) 14.5 33]);
end
figure('color','w')
hold all
scatter(kp2,logL5,100,'ko','markerfacecolor','w','linewidth',1)
plot([8.9 8.9],[0 logL5(kp2==8.9)],':','color',[.7 .7 .7],'linewidth',2)
text(8.9,-10000,'8.9','color',[.4 .4 .4])
[minlogL,I]=min(logL5);
plot([kp2(I) kp2(I)],[0 minlogL],'r:','linewidth',2)
text(kp2(I),-10000,num2str(kp2(I)),'color','r')
title('-Sum(log(likelihood)) nearby kp2(8.9)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kp3=14.5
kp3=0:14.5:700;
logL6=nan(numel(kp3),1);
for i=1:numel(kp3);
    [logL6(i)]=makeSSE5(upo,d,coh,pstd,[350 24 7.5 4.3 8.9 kp3(i) 33]);
end
figure('color','w')
hold all
scatter(kp3,logL6,100,'ko','markerfacecolor','w','linewidth',1)
plot([14.5 14.5],[0 logL6(kp3==14.5)],':','color',[.7 .7 .7],'linewidth',2)
text(14.5,-10000,'14.5','color',[.4 .4 .4])
[minlogL,I]=min(logL6);
plot([kp3(I) kp3(I)],[0 minlogL],'r:','linewidth',2)
text(kp3(I),-10000,num2str(kp3(I)),'color','r')
title('-Sum(log(likelihood)) nearby kp3(14.5)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc
%%
%kp4=33
kp4=0:33:700;
logL7=nan(numel(kp4),1);
for i=1:numel(kp4);
    [logL7(i)]=makeSSE5(upo,d,coh,pstd,[350 24 7.5 4.3 8.9 14.5 kp4(i)]);
end
figure('color','w')
hold all
scatter(kp4,logL7,100,'ko','markerfacecolor','w','linewidth',1)
plot([33 33],[0 logL7(kp4==33)],':','color',[.7 .7 .7],'linewidth',2)
text(33,-10000,'33','color',[.4 .4 .4])
[minlogL,I]=min(logL7);
plot([kp4(I) kp4(I)],[0 minlogL],'r:','linewidth',2)
text(kp4(I),-10000,num2str(kp4(I)),'color','r')
title('-Sum(log(likelihood)) nearby kp4(33)')
xlabel('K')
ylabel('-Sum(log(likelihood))')
box off
toc