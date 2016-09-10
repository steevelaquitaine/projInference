

%SLsimulateBayesWithBernoulliPrior(1,0.65,[45 135])

function SLsimulateBayesWithBernoulliPrior(kl,pMax,pDir)

figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [60 90 135 180 210];

%motion direction space
dir = 1:1:360;

%posterior for each evidence
for i = 1 : length(evidence)
    
    subplot(length(evidence),1,i)
    ylim([0 1])
    box off
    hold all
    
    %llh (circular)
    llh = vmPdfs(dir,evidence(i),kl,'norm');
    plot(dir,llh,'color',[.65 .65 .65],'linesmoothing','on')
    
    %prior
    prior = zeros(1,length(dir));
    prior(pDir) = [1-pMax pMax]; 
    prior=prior';
    plot(dir,prior,'b','linesmoothing','on')
    
    %posterior
    post = prior.*llh;
    post = post/sum(post);
    plot(dir,post,':r','linesmoothing','on','linewidth',3)
   
    %mark MAP estimates
    [~,mxp] = max(post);
    plot(mxp,1,'r.','markersize',10)
    
    %legend
    if i==1;
        legend('llh','prior','posterior','MAP')
        legend('boxoff')
        title('Bayesian estimation with Bernoulli(/circular) prior(/LLH)')
    end
end
