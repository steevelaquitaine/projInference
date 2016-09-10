

%SLsimulatePoissonPriorAndLLH.m
%
% author: steeve laquitaine
%   date: 141122
%purpose: simulate Bayesian posterior in an example trial when llh is a
%         Poisson distribution and prior is a von Mises.
%         
% usage:
%       SLsimulatePoissonPriorAndLLH
%
%Predictions:
% - estimates are attracted toward prior near prior and toward motion
%   distant from the prior.
% - But I can never get bimodality

function SLsimulatePoissonPriorAndLLH

figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [5 45 135 210];

%motion direction space
dir = 1:1:360;

%posterior for each evidence
maxY = [];
for i = 1 : length(evidence)
    
    ax(i) = subplot(length(evidence),1,i);
    box off
    hold all
    
    %llh
    llh = SLPoissonPdfs(dir,evidence(i),'norm');
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    prior = SLPoissonPdfs(dir,225,'norm');
    area(prior,'facecolor','b','linestyle','none')
    
    %posterior
    post = prior.*llh;
    post = post/sum(post);
    %plot(post,'--r','linesmoothing','on')
    area(post,'facecolor','r','linestyle','none')
    
    %legend
    if i==1;
        legend('llh','prior','posterior')
        legend('boxoff')
    end
    
    %axes
    maxY = [maxY;max([llh(:);prior(:);post(:)])];
    xlim([0 360])
end

%axes
set(ax,'ylim',[0 max(maxY)])

SLremoveDeadSpace

% Create textbox
annotation(gcf,'textbox',...
    [0.3 0.96 0.20 0.04],...
    'String','Posterior for 4 evidences produced by the same motion direction',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.3 0.8 0.2 0.04],...
    'String',['Trial motion evidence i = ',num2str(evidence(1))],...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.3 0.6 0.2 0.04],...
    'String',['Trial motion evidence i = ',num2str(evidence(2))],...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.3 0.4 0.2 0.04],...
    'String',['Trial motion evidence i = ',num2str(evidence(3))],...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.321663780663779 0.183563954729874 0.206142857142857 0.045238095238095],...
    'String',['Trial motion evidence i = ',num2str(evidence(4))],...
    'FitBoxToText','off',...
    'EdgeColor','none');

