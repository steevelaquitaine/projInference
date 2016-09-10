
%SLsimAveragePPCafterAdaptationSpikeRangeChange.m
%
%       $Id: SLsimAveragePPCafterAdaptationGainChange.m $
%        by: steeve laquitaine
%      date: 141220
%
%   purpose: Create repeated simulation of the outputs of a probabilistic
%            population code the spike range on which Poisson variability is defined changes:
%            - simulation of average population response
%            - simulation of average likelihood               
%           
%     usage:
% 
%          output = SLsimAveragePPCafterAdaptationSpikeRangeChange(1,90,...
%               4,1,360,225,10,0.7,[1:1:60],'display=off');
%
%
%Description: same as in SLsimPPCafterAdaptation.m
%
%
%References:
%
%    PPC
%       http://homepages.inf.ed.ac.uk/pseries/CCN14/lab4.pdf
%       Jazayeri et al,2006, Nat.Neu
%       Ma et al,2006, Nat.Neu
%
%    Adaptation
%       Benucci, A., Saleem, A. B. & Carandini, M. 
%       Adaptation maintains population homeostasis in primary 
%       visual cortex. Nature Publishing Group 16, 724?729 (2013).
%
function output = SLsimAveragePPCafterAdaptationSpikeRangeChange(numSim,stim,base,maxResp,...
    nbNeuron,Adap,kAdap,magAdap,rSpace,Optiondisplay)

t = tic;

%simulate for each adaptation gain
for i = 1 : size(rSpace,1)
    
    %this gain
    rSpace_i = rSpace(i,:);
    
    %clear
    clear o_i
    
    %simulate
    o_i = SLsimAveragePPCafterAdaptation(numSim,stim,base,maxResp,...
        nbNeuron,Adap,kAdap,magAdap,rSpace_i,Optiondisplay);    
    
    %backup
    output.diSpace                    = o_i.diSpace;
    output.D                          = o_i.D;
    output.meanLLHsim(:,i)            = o_i.meanLLHsim;
    output.meanPopRESPsim(:,i)        = o_i.meanPopRESPsim;
    output.modeMeanLLHsim(:,i)        = o_i.modeMeanLLHsim;
    output.Term_Jaz_meanSim(:,i)      = o_i.Term_Jaz_meanSim;
    output.Term_Pooledfi_meanSim(:,i) = o_i.Term_Pooledfi_meanSim;
    output.Term_Pooledni_meanSim(:,i) = o_i.Term_Pooledni_meanSim;
    output.LLHnotScaled(:,i)          = o_i.LLHnotScaled;
    output.PRgivenDi{i}               = o_i.PRgivenDi;
    
    %clear
    clear o_i    
end
 

%Plot correct (rSpace=0:1:60) vs incorrect Poisson noise (rSpace=1:1:61)
%-----------------------------------------------------------------------
%e.g., Poisson variability at 225 deg for neuron with preference for 225 deg
figure('color','w','position',[507 509 458 303])
exneur = 225;
exdir = 225;
barcol = SLlinspecer(10);
for i = 1 : size(rSpace,1) 
    hold all
    bar(rSpace(i,:),output.PRgivenDi{i}(:,exdir,exneur),...
        'facecolor',barcol(i,:),...
        'edgecolor','w')
    legend
    legi{i} = ['Spike range: ',num2str(i)];
    
    %stats
    meanResp(i) = SLmakeMean(rSpace(i,:),output.PRgivenDi{i}(:,exdir,exneur));
    stdResp(i) = SLmakeStd(rSpace(i,:),output.PRgivenDi{i}(:,exdir,exneur));
end
alpha(0.3)
axis tight
xlim([0 15])
box off
title({'Poisson noise for 225 deg-tuned neuron at 1 deg',['Mean. resp:[',...
    num2str(meanResp),']'],['Std. resp  :[',num2str(stdResp),']']})
legend(legi)
legend('boxoff')
xlabel('Spike count (spikes)')
ylabel('Probability')


%plot (log form of each of the three terms)
%------------------------------------------
%Likelihood formulation based on Jazayeri et al., 2006 Nature N.
%term1
figure('color','w','position', [456 270 1060 297])
h(1) = subplot(1,5,1);
plot(output.diSpace,output.Term_Jaz_meanSim)
box off
axis tight
title({'Term 1: LLH approx (J. et al)','$\exp{\sum_{i=1}^N (n_i \times \log{f_i(\theta)}))}$'},'interpreter','latex')
ylabel('Log10 scaled (au)')

%term2
h(2) = subplot(1,5,2);
plot(output.diSpace,output.Term_Pooledfi_meanSim)
box off
axis tight
title({'Term 2',': $\exp{(- \sum_{i=1}^N (f_i(\theta)))}$'},'interpreter','latex')

%term3
h(3) = subplot(1,5,3);
plot(output.diSpace,output.Term_Pooledni_meanSim)
box off
axis tight
minP = min(min(output.Term_Pooledni_meanSim));
maxP = max(max(output.Term_Pooledni_meanSim));
ylim([0.9*minP 1.1*maxP])
title({'Term 3:','$\exp{(- \sum_{i=1}^N (\log(n_i!)))}$'},'interpreter','latex')
xlabel('Motion. direction (degree)')

%LLH (not scaled)
h(4) = subplot(1,5,4);
semilogy(output.diSpace,output.LLHnotScaled)
box off
axis tight
title({'LLH (not scaled)','term1*term2*term3'},'interpreter','latex')

%LLH (scaled to 1)
h(5) = subplot(1,5,5);
Z = sum(output.LLHnotScaled);
Zall = Z(ones(output.D,1),:);
output.LLHScaled = output.LLHnotScaled./Zall;
plot(output.diSpace,output.LLHScaled)
box off
axis tight
title('LLH (scaled)','interpreter','latex')
ylabel('Proba)')


%legend
for i = 1 : size(rSpace,1)
    leg_i{i} = ['Spike range: ',num2str(i)];
end
legend(leg_i)
legend('boxoff')


% SLConventionUp
SLremoveDeadSpace

toc(t)


