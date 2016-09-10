

%SLsimulateBayesVsCompBimodalPrior.m
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
%
%       output = SLsimulateBayesVsCompBimodalPrior(2.98,[165 285],3.71,'experimentalPrior')
%
%Predictions:
%
%Difference in models predictions is maximized when:
%- motion is weak enough to exploit prior
%- prior is weak enough to produce direction distant from the modes but not
%strong enough to compete with motion.


function output = SLsimulateBayesVsCompBimodalPrior(kllh,modesPrior,kprior,varargin)

%call for help
if ieNotDefined('kllh')
    help SLsimulateBayesVsCompBimodalPrior
    return
end

figure('color','w')
clf

%motion direction of an actual experiment are chosen
%Certain directions from the sample direction won't be 
%displayed on screen. We only look at the motion directions
%that will actually be displayed.
%prior (for actual experiment)
if strcmp(varargin{:},'experimentalPrior')
    o = SLmakeDiscreteMixtureVM(kprior,5:10:355,modesPrior,107);
    evidence = o.parameter.dir.sample.degree(o.parameter.dir.count~=0);
else 
    evidence = [5:10:355];
end

%motion direction space
dir = 1:1:360;

%posterior for each evidence
for i = 1 : length(evidence)
    
    subp(i) = subplot(length(evidence),1,i);
    ylim([0 0.06])
    box off
    hold all
    
    %llh
    output.llh = vmPdfs(dir,evidence(i),kllh,'norm');
    area(output.llh,'facecolor',[.65 .65 .65],'linestyle','none')
    
    %prior
    output.prior = 0.5*vmPdfs(dir,modesPrior(1),kprior,'norm') + 0.5*vmPdfs(dir,modesPrior(2),kprior,'norm');
    area(output.prior,'facecolor','b','linestyle','none')
    
    %posterior
    post = output.prior.*output.llh;
    output.post = post/sum(post);
    area(output.post,'facecolor','r','linestyle','none')
    
    %std posterior
    output.stdpost(i) = SLmakeStd(dir,post);
    
    %graphics
    alpha(0.7)
    
    %readouts
    %switching
    maxSwitc = unique(max(max([output.prior output.llh])));
    hold all
    plot(modesPrior(1),maxSwitc,'ob','markersize',6,'markerfacecolor','b')
    plot(modesPrior(2),maxSwitc,'ob','markersize',6,'markerfacecolor','b')
    plot(evidence(i),maxSwitc,'ob','markersize',6,'markerfacecolor','b')
    
    %Bayes MAP
    output.postApprox = ceil(output.post*10^7)/10^7;
    maxpost = max(output.postApprox);
    argmax = find(output.postApprox==maxpost);
    for j = 1 : numel(argmax)
        plot(argmax(j),maxpost,'or','markersize',6,'markerfacecolor','r')
    end
    
    %legend
    if i==1;
        legend('llh','prior','posterior','Switching readout')
        legend('boxoff')
    end
    ylim([0 max([output.llh(:); output.prior(:); output.post(:)])])
    xlim([0 360])
end

%Labels
for i = 1 : length(evidence)
    ylabel(subp(i),'Probability')
end
xlabel('Motion directions')
% SLConventionUp
