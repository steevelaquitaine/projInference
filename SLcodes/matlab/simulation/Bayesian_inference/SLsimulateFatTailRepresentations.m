

%SLsimulateFatTailRepresentations.m
%
% author: steeve laquitaine
%   date: 141122
%purpose: simulate Bayesian posterior in an example trial when prior is a
%         von Mises with a fat tail (mixture of von Mises and uniform)
%
%reference: Stocker and Simoncelli, 2006, NN

function SLsimulateFatTailRepresentations


figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [90 135 180 210];

%motion direction space
dir = 1:1:360;

%posterior for each evidence
for i = 1 : length(evidence)
    
    subplot(length(evidence),1,i)
    ylim([0 0.06])
    box off
    hold all
    
    %llh
    llh = vmPdfs(dir,evidence(i),80,'norm');
    weightTail = 0;
    uniform = (1/360)*ones(numel(dir),1);
    llh = (1 - weightTail)*llh  + weightTail*uniform;
    plot(llh,'color',[.65 .65 .65],'linesmoothing','on')
    
    %prior
    prior = vmPdfs(dir,225,40,'norm');
    weightTail = 0.00001;
    prior = (1 - weightTail)*prior  + weightTail*uniform;
    plot(prior,'b','linesmoothing','on')
    
    %posterior
    post = prior.*llh;
    post = post/sum(post);
    plot(post,'--r','linesmoothing','on')
   
    %legend
    if i==1;
        legend('llh','prior','posterior')
        legend('boxoff')
    end
end
SLremoveDeadSpace