
%bugs
%why does -logL not vary with kl ???


tic
%Look at how log likelihood of any given data behaves as we change the
%parameters: data, displayed directions and width of the prior.
%658 seconds
data=1:50:360;
displ=5:10:355; 
kp=0.001:50:699;
kl=0.001:1:699;
for a=1:numel(data)
    for b=1:numel(displ)
        for c=1:numel(kp)
            for d=1:numel(kl)
                logL(a,b,c,d)=GirshickML(data(a),displ(b),kp(c),kl(d));
                minuSumlogL(c,d)=-sum(sum(logL(:,:,c,d))); 
            end
        end
    end
    toc
end            

%% plot the likelihood of entire data set for each parameter combination
figure('color','w')
subplot(2,3,1)
imagesc(kp,kl,minuSumlogL)
set(gca,'YDir','normal')
xlabel('klikelihood')
ylabel('kprior')
colorbar
title('Colder colors indicate max likelihood')
box off

%Examples of the influence of K prior on the -log likelihood
subplot(2,3,5)
hold all
plot(kp,minuSumlogL(:,kl==7.001),'-o')
xlabel('Kprior')
ylabel('-Log likelihood')
ylim([1.5 3*10^4])


%Examples of the influence of K prior on the -log likelihood
subplot(2,3,6)
hold all
plot(kl,minuSumlogL(7,:),'-o')
xlabel('Klikelihood')
ylabel('-Log likelihood')
ylim([1.5 3*10^4])


%data is far from displayed
logL1=-GirshickML(1,180,1,7);
logL2=-GirshickML(1,180,300,7);
logL3=-GirshickML(1,180,699,7);
plot([1 300 699],[logL1 logL2 logL3],'-o')

%data is closed to displayed
logL1=-GirshickML(1,100,1,7);
logL2=-GirshickML(1,100,300,7);
logL3=-GirshickML(1,100,699,7);
plot([1 300 699],[logL1 logL2 logL3],'--o')

xlabel('Kprior')
ylabel('-Log likelihood')
xlim([0 700])
ylim([1.5 25])
box off
title('kp does not influence -logL much')


%Examples of the influence of K likelihood on the -log likelihood
subplot(2,3,3)
hold all

%data is far from displayed
logL1=-GirshickML(1,180,100,1);
logL2=-GirshickML(1,180,100,300);
logL2=-GirshickML(1,180,100,300);
logL3=-GirshickML(1,180,100,699);
plot([1 300 699],[logL1 logL2 logL3],'-o')

%data is closed to displayed
logL1=-GirshickML(1,1,100,1);
logL2=-GirshickML(1,1,100,300);
logL3=-GirshickML(1,1,100,699);
plot([1 300 699],[logL1 logL2 logL3],'--o')
xlabel('Klikelihood')
ylabel('-Log likelihood')
xlim([0 700])
ylim([1.5 25])
box off
title('kl influences strongly -logL')
h=legend('disp:180 - estimated:1','disp:1 - estimated:1')
set(h,'box','off','location','Best')
