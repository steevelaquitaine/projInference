

%SLsimulateFatTailPrior.m
%
% author: steeve laquitaine
%   date: 141122
%purpose: simulate Bayesian posterior in an example trial when prior is a
%         von Mises with a fat tail (mixture of von Mises and uniform)
%
%reference: Stocker and Simoncelli, 2006, NN
%
%
%usage:
%       SLsimulateFatTailPrior(0.8,2.9,0.7)
%       SLsimulateFatTailPrior(0.1,4,0.7)
%       SLsimulateFatTailPrior(0.6,15,33)
%       output = SLsimulateFatTailPrior(.6,2.98,33)
%
%Predictions:
% - This model often produces posteriors with two peaks.
% - It most often predicts bimodal estimates distribution with one peak in
%   between prior and llh and the other at the likelihood.
% - It matches the data best when prior is very strong (600) compared to
%   likelihood (1) and have very fat tails (weightTail=0.95)
% - It may be able to explain the data when readout is sampling but
%   predictions are complex. They change with llh distance to the prior and
%   bimodality requires very fat tail strong priors. Intuitively it would
%   indicates very poor learning for directions distant from the prior with
%   maximal learning at the prior.
%
%   SLsimulateFatTailPrior(0.95,1,600)



function output = SLsimulateFatTailPrior(weightTail,kl,kp)

%call for help
if ieNotDefined('weightTail')
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
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    prior = vmPdfs(dir,225,kp,'norm');
    uniform = (1/360)*ones(numel(dir),1);
    prior = (1 - weightTail)*prior  + weightTail*uniform;
    area(prior,'facecolor','b','linestyle','none')
    
    %posterior
    post = prior.*llh;
    post = post/sum(post);
    area(post,'facecolor','r','linestyle','none')
    
    %std posterior
    output.stdpost(i) = SLmakeStd(dir,post);
    
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

