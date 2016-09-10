
%SLsimAveragePPCafterAdaptationGainChange.m
%
%       $Id: SLsimAveragePPCafterAdaptationGainChange.m $
%        by: steeve laquitaine
%      date: 141216
%
%   purpose: Create repeated simulation of the outputs of a probabilitic
%            population code when adaptation gain changes:
%            - simulation of average population response
%            - simulation of average likelihood               
%           
%     usage:
% 
%          output = SLsimAveragePPCafterAdaptationGainChange(10,90,2,1,360,225,10,[0 0.6 0.7],'display=off');
%
%Description: same as in SLsimPPCafterAdaptation.m
%
%
%
%Predictions:    
%
%   - In the absence of adaptation, LLH peak at the adapter > LLH peak at 
%     the motion direction.
%
%   - When adaptation is strong, LLH peak at the adapter > LLH peak at 
%     the motion direction.
%
%     output = SLsimAveragePPCafterAdaptationGainChange(10,90,1,1,360,225,10,...
%          [0 0.5 0.7],'display=off');
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
function output = SLsimAveragePPCafterAdaptationGainChange(numSim,stim,base,maxResp,...
    nbNeuron,Adap,kAdap,magAdap,Optiondisplay)

t = tic;

%simulate for each adaptation gain
for i = 1 : numel(magAdap)
    
    %this gain
    magAdap_i = magAdap(i);
    
    %clear
    clear o_i
    
    %simulate
    o_i = SLsimAveragePPCafterAdaptation(numSim,stim,base,maxResp,...
        nbNeuron,Adap,kAdap,magAdap_i,Optiondisplay);    
    
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
    
    %clear
    clear o_i    
end

%plot (LLH three terms are in their log form for visibility)
%-----------------------------------------------------------
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
for i = 1 : numel(magAdap)
    leg_i{i} = ['Adap. Gain: ',num2str(magAdap(i))];
end
legend(leg_i)
legend('boxoff')


% SLConventionUp
SLremoveDeadSpace

toc(t)


