%test

% figure('color','w');
clf
subplot(3,1,1)
hold all
motion=1:1:360;

%Encoding of prior and likelihood
%................................
%llhs
llh80=vmPdfs(motion,90,80,'norm');
plot(motion,llh80,'b--','linesmoothing','on')
text(90,1.3*max(llh80),'llhs','fontsize',14)

%priors
prior80=vmPdfs(motion,225,80,'norm');
prior40=vmPdfs(motion,225,40,'norm');
plot(motion,prior80,'r','linesmoothing','on')
plot(motion,prior40,'color',[0.1 0 0],'linesmoothing','on')
text(225,1.3*max(prior80),'Priors','fontsize',14)
xlim([0 360])
ylabel({'Probabilistic','representations'},'fontsize',14)
box off


%Competition between prior and llh
%.................................
subplot(3,1,2)
hold all
%calculate P(llh)=softmax(stdllh-stdPrior)
%Does "stdPrior-stdllh" equals std"Prior-llh"?

%Selection by mutual inhibition between llh and prior representations
%We use log(llh) as in (Jazayeri et al,2006) - log(prior) 
%log(llh)-log(Prior) looks like predictive coding's prediction error..
%prior representations inhibit llh representations.
%note: integration would be log(llh) + log(prior)

%strong prior 
compet=(llh80./prior80)/sum(llh80./prior80);
plot(motion,compet,'--','color',[1 0 1],'linesmoothing','on')
text(90,max(compet),'llh/prior','fontsize',14)

compet=(prior80./llh80)/sum(prior80./llh80);
plot(motion,compet,'color',[1 0 1],'linesmoothing','on')
text(225,max(compet),'prior/llh','fontsize',14)

%weak prior
compet=(llh80./prior40)/sum(llh80./prior40);
plot(motion,compet,'--','color',[.1 0 1],'linesmoothing','on')

compet=(prior40./llh80)/sum(prior40./llh80);
plot(motion,compet,'color',[.1 0 1],'linesmoothing','on')

box off
xlim([0 360])
ylabel('Competition','fontsize',14)
xlabel('Motion direction (degrees)','fontsize',14)
set(gcf,'renderer','painters');


%then the resulting representation is pooled 
%(integration over motion direction units)
%...........................................
subplot(3,1,3)




