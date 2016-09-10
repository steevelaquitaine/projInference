
%flatTailVonMisesPrior
%We obtain bimodal posterior only if both llh and prior have flat tails.
%This means that subjects overweight the probability that direction further
%away from the motion direction and the most likely direction occur.

figure('color','w')
subplot(4,2,1)
ullh=135;
unic=(1./360)*ones(360,1);

%Bayesian inference
llh=vmPdfs(1:1:360,ullh,2.6,'norm');
prior=vmPdfs(1:1:360,225,0.8,'norm');
posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)
title('Bayesian inference with Von mises','fontsize',14)
title('Bayesian inference with fat tail von mises','fontsize',14,...
    'fontweight','bold')
ylabel('Probability','fontsize',14)


%Bayesian inference
subplot(4,2,3)
llh=vmPdfs(1:1:360,ullh,2.6,'norm');
prior=vmPdfs(1:1:360,225,2.7,'norm');
posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);
hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)



%Bayesian inference
subplot(4,2,5)
llh=vmPdfs(1:1:360,ullh,2.6,'norm');
prior=vmPdfs(1:1:360,225,14.5,'norm');
posterior=llh.*prior;
posterior=posterior/sum(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)


%Bayesian inference
subplot(4,2,7)
llh=vmPdfs(1:1:360,ullh,2.6,'norm');
prior=vmPdfs(1:1:360,225,160,'norm');
posterior=llh.*prior;
posterior=posterior/sum(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)
xlabel('Motion direction','fontsize',14)


%Bayesian inference with fat tail von mises
subplot(4,2,2)
llhc=vmPdfs(1:1:360,ullh,2.6,'norm');
llh=0.5.*(llhc+unic);

priorc=vmPdfs(1:1:360,225,0.8,'norm');
prior=0.5.*(priorc+unic);

posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)
title('Bayesian inference with Von mises','fontsize',14)
title('Bayesian inference with fat tail von mises','fontsize',14,...
    'fontweight','bold')



%Bayesian inference with fat tail von mises
subplot(4,2,4)
llhc=vmPdfs(1:1:360,135,2.6,'norm');
llh=0.5.*(llhc+unic);

priorc=vmPdfs(1:1:360,225,2.7,'norm');
prior=0.5.*(priorc+unic);

posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(135,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)



%Bayesian inference with fat tail von mises
subplot(4,2,6)
llhc=vmPdfs(1:1:360,ullh,2.6,'norm');
llh=0.5.*(llhc+unic);

priorc=vmPdfs(1:1:360,225,14.5,'norm');
prior=0.5.*(priorc+unic);

posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)


%Bayesian inference with fat tail von mises
subplot(4,2,8)
llhc=vmPdfs(1:1:360,ullh,2.6,'norm');
llh=0.5.*(llhc+unic);

priorc=vmPdfs(1:1:360,225,160,'norm');
prior=0.5.*(priorc+unic);

posterior=llh.*prior;
posterior=posterior/sum(posterior);
[maxp,MAP]=max(posterior);

hold all
plot(llh,'--','linesmoothing','on')
text(ullh,1.3*max(llh),'llh','color','b','fontsize',14)
plot(prior,'k','linesmoothing','on')
text(225,1.3*max(prior),'prior','color','k','fontsize',14)
plot(posterior,'r','linesmoothing','on')
text(MAP,1.3*max(posterior),'posterior','color','r','fontsize',14)
xlabel('Motion direction','fontsize',14)

