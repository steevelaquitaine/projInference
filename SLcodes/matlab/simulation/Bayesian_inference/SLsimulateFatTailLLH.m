

%SLsimulateFatTailLLH.m
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
%       SLsimulateFatTailLLH(0.01,7,0.7)
%       SLsimulateFatTailLLH(0.1,4,0.7)
%       SLsimulateFatTailLLH(0.3,4,2)
%       SLsimulateFatTailLLH(0.1,4,33)
%       output = SLsimulateFatTailLLH(.2,4,10)


function output = SLsimulateFatTailLLH(weightTail,kl,kp)

%call for help
if ieNotDefined('weightTail')
    help SLsimulateFatTailPrior
    return
end

figure('color','w')
clf

%four sensory evidences are sampled in four trials
evidence = [55 135 180 210];

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
    uniform = (1/360)*ones(numel(dir),1);
    llh = (1 - weightTail)*llh  + weightTail*uniform;
    area(llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    prior = vmPdfs(dir,225,kp,'norm');
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

output.llh = llh;

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

