

%SLsimulatePTweightedPriorAndLLH.m
%
% author: steeve laquitaine
%   date: 141122
%purpose: simulate Bayesian posterior in an example trial when prior is a
%         von Mises weighted by Prospect theory weighting function.
%
%reference: 
%           Tversky and kahneman, 1981,
%           Tversky and Fox 1995,
%           Neuroeconomics book, page 158
%
%
%usage:
%       SLsimulatePTweightedPriorAndLLH(1,.1,2,4.77)
%
%Predictions:
%
% similar to von mises llh and prior. Unimodal posterior.
% cannot explain the bimodal data observed when prior and llh are weak 
% because tails are not heavy enough.



function SLsimulatePTweightedPriorAndLLH(delta,gamma,kl,kp)

%call for help
if ieNotDefined('delta')
    help SLsimulateFatTailPrior
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
    llh = delta.*llh.^gamma./(delta.*llh.^gamma + (1 - llh).^gamma);
    llh = llh./sum(llh);
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    %exact prior
    prior_exact = vmPdfs(dir,225,kp,'norm');

    %Prosp. Theory weighted prior (Tversky and Fox, 1995)
    prior = delta.*prior_exact.^gamma./(delta.*prior_exact.^gamma + (1 - prior_exact).^gamma);
    prior = prior./sum(prior);
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


%plot
%compare prospect theory weighted prior with exact prior
figure;
hold all;
area(prior_exact,'facecolor',[.5 .5 .5],'linestyle','none')
area(prior,'facecolor','b','linestyle','none')
alpha(0.3)
legend('Exact prior','PT weighted prior')
legend show

