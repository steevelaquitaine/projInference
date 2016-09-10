

%SLsimulateLaplacePrior.m
%
% author: steeve laquitaine
%   date: 141122
%purpose: simulate Bayesian posterior in an example trial when prior is a
%         Laplace distribution
%   note: this is not a circular distribution but with reasonable
%   parameters we can get valid insight about the shape of the posterior.
%
%reference: Stocker and Simoncelli, 2006, NN
%
%
%usage:
%       SLsimulateLaplacePrior(3,100)
%
%Predictions:
% - This model always produces posterior with one peak.



function SLsimulateLaplacePrior(kl,b)

%call for help
if ieNotDefined('kl')
    help SLsimulateLaplacePrior
    return
end

figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [5 135 180 210];

%motion direction space
dir = 1:1:360;

%posterior for each evidence
for i = 1 : length(evidence)
    
    subplot(length(evidence),1,i)
    ylim([0 0.06])
    box off
    hold all
    
    %llh
    llh = vmPdfs(dir,evidence(i),kl,'norm');
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    prior = SLLaplacePdf(dir,225,b,'norm');
    area(prior,'facecolor','b','linestyle','none')
    
    %posterior
    post = prior.*llh;
    post = post/sum(post);
    area(post,'facecolor','r','linestyle','none')
    
    %legend
    if i==1;
        legend('llh','prior','posterior')
        legend('boxoff')
    end
    ylim([0 max([llh(:);prior(:);post(:)])])
    xlim([0 360])
end

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

