
function getLogL
%bugs
%why does -logL not vary with kl ???
tic
load data

%Look at how log likelihood of any given data behaves as we change the
%parameters: data, dayed directions and width of the prior.
%658 seconds
%data=1:50:360;
%d=5:10:355;
kp=1:1:30;

%kl obtained from fitting
kl=2;

for i=1:numel(kp)
    for j=1:numel(upo)
        logL(j,i)=GirshickML(upo(j),d(j),kp(i),kl);
        minuSumlogL(i)=-sum(logL(:,i));
    end
    hold all
    plot(i,minuSumlogL(i),'o')
    drawnow
    toc
end
keyboard
