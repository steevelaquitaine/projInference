

%SLsimulateBayesFreePrior_CircLLH(0.5,[45 90 135],[0.25 0.5 0.25])

function SLsimulateBayesFreePrior_CircLLH(kl,pDir,pP)

figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [20 60 110 135 180];

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
    prior(pDir) = pP; 
    
    %warning
    if sum(pP)~=1
        disp('Prior probabilities pP must sum to 1');
        return
    end
    
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
        title('Bayesian estimation free(/circular) prior(/LLH)')
    end
end