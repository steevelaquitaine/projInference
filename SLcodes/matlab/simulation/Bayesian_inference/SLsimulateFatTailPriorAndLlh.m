

%SLsimulateFatTailPriorAndLlh.m
%
%   author: steeve laquitaine
%     date: 141122
%  purpose: simulate Bayesian posterior in an example trial when prior and llh are 
%           von Mises with a fat tail (mixture of von Mises and uniform)
%
%reference: Stocker and Simoncelli, 2006, NN
%
% usage:
%       SLsimulateFatTailPriorAndLlh(0.001,5.4 , 5.4)
%       SLsimulateFatTailPriorAndLlh(0.3, 10^12, 15)
%       SLsimulateFatTailPriorAndLlh(0.9, 24.5, 5.4)
%       SLsimulateFatTailPriorAndLlh(0.5,4,4400)   
%
%Predictions:
%
% - This model often produces a posterior with two peaks.
% - It predicts unimodal estimates distributions for most of the parameters
% - When Readout is MAP it predicts bimodal estimates distributions only 
%   when prior and llh strengths are equal (peaks have equal magnitudes)
% - When Readout is Sampling it predicts bimodal estimates distributions
%   and the peaks at the prior and llh depends on prior and llh strength as
%   observed in the data. When tails are fatter, the model predicts more
%   random estimation of directions that are not llh or prior mean (lapse 
%   rate - like). This is not seen in the data. So Tails must be very thin.
% 
% - Intuitively it would indicate that sensory inputs are corrupted 
%   by non Gaussian noise that may corrupt the learning of the prior or
%   that all representations are skewed toward the most likely signal.


function SLsimulateFatTailPriorAndLlh(weightTail,kl,kp)


%%
figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [5 20 40 155 180 210];

%motion direction space
dir = 1:1:360;

%posterior for each evidence
for i = 1 : length(evidence)
    
    subplot(length(evidence),1,i)
    box off
    hold all
    
    %llh
    llh = vmPdfs(dir,evidence(i),kl,'norm');
    uniform = (1/360)*ones(numel(dir),1);
    llh = (1 - weightTail)*llh  + weightTail*uniform;
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    prior = vmPdfs(dir,225,kp,'norm');
    prior = (1 - weightTail)*prior  + weightTail*uniform;
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


