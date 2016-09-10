

  %Author: steeve laquitaine
    %date: 140424
   %usage: AdaptiveGainBayesianLearning
 %purpose: Simulate change of posterior belief about motion direction
           %during learning.
           %Two models of learning are used:
           % - Bayesian inference (ideal optimal observer)
           % - Adaptive gain control (Summerfield et al., 2013, Neuron)

%Question: Does Adaptive gain control explains why subjects categorize
           %probabilistically between the motion direction (llh mean) and 
           %the most likely direction (prior mean)?
           %My guess is that adaptive gain rescale the gaussian
           %distribution in stronger non-gaussian distribution that are
           %integrated by Bayesian operations. Behavior thus deviates from
           %optimality because adaptive gain distorts llh representation
           %which distorts the learnt prior. When the two are combined by
           %bayesian operation they do not produce a posterior estimate
           %intermediate to prior and likelihood but switch between llh
           %mean and prior mean. 

function AdaptiveGainBayesianLearning

%%

%ideal Bayesian learning
Posterior(:,1)=vmPdfs(1:1:360,225,0,'norm');

%prior over motion directions
numTrials=1000;
Prior=vmPdfs(1:1:360,225,2.75,'norm');
motiondir=randsample(1:1:360,numTrials,true,Prior);

%llh
kllh=2.75;
llh=vmPdfs(1:1:360,motiondir,ones(numTrials,1)*kllh,'norm');

%Update posterior (learning)
for t=2:numTrials;
  
    %Ideal observer as a control of optimality
    %update posterior at time t with new scaled llh
    Posterior(:,t)=Posterior(:,t-1).*llh(:,t);
    
    %scale to probabilities (sum to 1)
    Posterior(:,t)=Posterior(:,t)./sum(Posterior(:,t));
end

%plot posteriors
figure('color','w')
title(['Prior (k=',num2str(kllh),'); ',...
    num2str(numTrials),' trials'])
ylabel('Probability')
xlabel('Motion directions (degrees)')
for t=1:1000; 
    plot(Posterior(:,t),'linesmoothing','on','linewidth',2)
    ylim([0 max(Posterior(:))])
    xlim([0 360])
    box off 
    drawnow
end


%With adaptive gain theory
%ScaledlogPosterior(:,1)=vmPdfs(1:1:360,225,0,'norm');
% for t=
    %scaling parameters of sigmoid function (Summerfield et al., 2013)
%     xk=
    
    %llh
%     llh(:,t)=(1+exp(-llh(:,t))).^-1;

   %Adaptive gain inference
    %scale llh by a sigmoid tranfer function (Summerfield et al., 2013)
%     scaledllh(:,t)=(1+exp(-llh(:,t))).^-1;
    
    %With adaptive gain. Each new llh is scaled
%     ScaledlogPosterior(:,t)=ScaledlogPosterior(:,t-1) + log(scaledllh(:,t)) ;
% end

