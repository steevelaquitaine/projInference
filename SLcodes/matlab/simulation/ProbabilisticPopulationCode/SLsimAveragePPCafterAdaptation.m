

%SLsimAveragePPCafterAdaptation.m
%
%       $Id: SLsimAveragePPCafterAdaptation.m$
%        by: steeve laquitaine
%      date: 141108
%
%   purpose: Create repeated simulations of the outputs of a probabilistic
%            population code:
%            - simulation of average population response
%            - simulation of average likelihood                
%           
%     usage:
%
%          output = SLsimAveragePPCafterAdaptation(10,90,1,1,360,225,4.7,1,[0:1:60],'display=on');
% 
%
%Description: see SLsimPPCafterAdaptation
%
%
%Predictions:    
%
%     Strong motion (maxResp = 60 spikes)
%       output = SLsimPPCafterAdaptation(180,60,1,100,225,4.7,1,[1:1:61],'display=on')
%       - LLH peaks at input direction
%       - LLH is strong (low std)
%
%     Weak motion (maxResp = 2 spikes)
%       output = SLsimPPCafterAdaptation(180,2,1,100,225,4.7,1,[1:1:61],'display=on')
%       - LLH shows two peaks: at input direction and at adapter
%       - the strongest the adaptation amplitude magAdap and the strongest
%       the adaptation.
%
%      It is because of the second term in Jazayeri's formulation of the
%      likelihood (sum of the tuning curves for each direction)
%
%
%
%Understanding why we get a peak at the adapter.
%Getting a peak at the adapter means that LLH slope increases
%as we approach the adapter (moving from 0º toward the adapter).
%LLH slope near the adapter depends on the slope of the three term.
%
%The slope of term3 is always constant across motion direction. So this
%term is irrelevant to LLH shape.
%
%The slope of Jazayeri's term actually decreases when approaching the
%adapter because tuning curves amplitude decreases when approaching the 
%adapter which pulls down Jazayeri's term.
%
%So term2: the pooled fi has to have a larger slope than Jazayeri's term 
%to create the positive slope observed when approaching the adapter.
%This is what we observe.
%
%A = [output.Term_Pooledfi_meanSim(200) output.Term_Pooledfi_meanSim(202) output.Term_Pooledfi_meanSim(225)];
%B = [output.Term_Jaz_meanSim(200) output.Term_Jaz_meanSim(202) output.Term_Jaz_meanSim(225)];
%
%Those are the values produced by averaging responses over 10 simulations
%A = [0.0021 0.0045 0.7453] * 1e-133
%B = [0.6159 0.3537 0.0508] * 1e-40
%
%to be continued
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
% note: works fine

function output = SLsimAveragePPCafterAdaptation(numSim,stim,base,maxResp,...
    nbNeuron,Adap,kAdap,magAdap,rSpace,Optiondisplay)

%run n simulations
output.LLH_Jazayeri = nan(360,numSim);
LLHsim = nan(360,numSim);
PopRESPsim = nan(nbNeuron,numSim);
for i = 1 : numSim
    
    %status
    sprintf(['sim ',num2str(i)])
    
    %simulation
    output = SLsimPPCafterAdaptation(stim,base,maxResp,nbNeuron,Adap,kAdap,magAdap,...
        rSpace,'display=on');
    
    %store
    LLHsim(:,i) = output.LLH_Jazayeri;
    PopRESPsim(:,i) = output.PopResp;
end

%store
output.LLHsim = LLHsim;
output.PopRESPsim = PopRESPsim;

%average LLH
output.meanLLHsim = mean(output.LLHsim,2);
output.meanPopRESPsim = mean(output.PopRESPsim,2);

%readout LLH (mode)
[maxLLH,mode] = max(output.meanLLHsim);
output.modeMeanLLHsim = output.diSpace(mode);


%-------------------------
%Plot steps of the process
%-------------------------
%Reasoning
%---------
%logLLH = sum_i( n_i * log(f_i(theta) ) - sum_i( f_i(theta) ) - sum_i(log(n_i!)));
%
%For clarity, let's define:
%
%A = sum_i( n_i * log(f_i(theta) )
%B = sum_i( log(f_i(theta)) )
%C = sum_i( log(n_i!)) )
%
%such that:
%
%logLLH = A - B - C
%
%Now,
%LLH not scaled = exp(logLLH)
%               = exp(A - B - C)
%               = exp(A) * exp(-B) * exp(-C)
%
%note: scaling does not change the shape of LLH which is 
%what we want to explain, so for clarity we do not scale yet.
%
%To see how each term contribute to LLH, we just plot 
%the separetely.
%
%note: exp(A) is Jazayeri approximated LLH so i'll call it 
%Term_Jaz =  exp(A)
%Term_Pooledfi = exp(-B)
%Term_Pooledni = exp(-C)
%
%
%Plot Jazayeri approx. of LLH
%
% A = sum_i( n_i * log(f_i(theta) );
% Term_Jaz =  exp(A);
%
%term1: Jazayeri approx. of LLH
n_i = output.meanPopRESPsim;
n_theta_i = repmat(n_i',output.D,1);
f_theta_i = output.f;

%This is jazayeri's formulation. For now f_theta_i cannot be too high ()be
%careful base and max response change because of numerical limits when 
%exponentiating large number.Numerical limits depends on the number of neurons
%e.g., term = + inf when {base = 26, nbNeurons = 10}
%                or when {base =  2, nbNeurons = 360}
%when we run output = SLsimAveragePPCafterAdaptationBaseChange(1,90,2,1,360,225,10,0.7,'display=off');
% A = sum(n_theta_i.*log(f_theta_i),2);
% Term_Jaz =  exp(A);
% 
% %term2 pooled tuning responses
% B = sum(f_theta_i,2);
% Term_Pooledfi = exp(-B);
% 
% %term3 pooled stim responses.This term = 0 when responses are high 
% %(e.g., base =4)if we increase base response and max response change.But it
% %doesn't matter as this term does not affect the shape of the LLH.
% C = sum(log(gamma(n_theta_i+1)),2);
% Term_Pooledni = exp(-C);
% 
% %LLH
% LLHnotScaled = Term_Jaz .* Term_Pooledfi.* Term_Pooledni;
% 
% %backup
% output.Term_Jaz_meanSim = Term_Jaz;
% output.Term_Pooledfi_meanSim = Term_Pooledfi;
% output.Term_Pooledni_meanSim = Term_Pooledni;
% output.LLHnotScaled = LLHnotScaled;
% 
% %plot
% if strcmp(Optiondisplay,'display=on')
%     
%     figure('color','w','position',[224 276 1060 297])
%     h(1) = subplot(1,6,1);
%     plot(Term_Jaz)
%     box off
%     axis tight
%     title({'Term 1: Jaz.approx of LLH','$\exp{\sum_{i=1}^N (n_i \times \log{f_i(\theta)}))}$'},'interpreter','latex')
%     
%     h(2) = subplot(1,6,[2 3 4]);
%     plot(Term_Pooledfi)
%     box off
%     axis tight
%     title('Term 2: $\exp{(- \sum_{i=1}^N (f_i(\theta)))}$','interpreter','latex')
%     xlabel('Motion. direction (degree)')
%     
%     h(3) = subplot(1,6,5);
%     plot(Term_Pooledni)
%     box off
%     axis tight
%     ylim([0.9*min(Term_Pooledni) 1.1*max(Term_Pooledni)])
%     title('Term 3: $\exp{(- \sum_{i=1}^N (\log(n_i!)))}$','interpreter','latex')
%     
%     h(4) = subplot(1,6,6);
%     plot(LLHnotScaled)
%     box off
%     axis tight
%     title('LLH (not scaled) term1*term2*term3','interpreter','latex')
%     
%     SLremoveDeadSpace
% end
% 

%This is jazayeri's formulation. I keep the term in their log form to
%prevent numerical limits with exponential when f_theta_i 
%base and max responses are large because of numerical limits when 
%exponentiating large number. Numerical limits of exponential also 
%depends on the number of neurons
%e.g., term = + inf when {base = 26, nbNeurons = 10}
%                or when {base =  2, nbNeurons = 360}
%when we run output = SLsimAveragePPCafterAdaptationBaseChange(1,90,2,1,360,225,10,0.7,'display=off');
Term_Jaz = sum(n_theta_i.*log(f_theta_i),2);

%term2 pooled tuning responses
Term_Pooledfi = sum(f_theta_i,2);

%term3 pooled stim responses.This term = 0 when responses are high 
%(e.g., base =4)if we increase base response and max response change.But it
%doesn't matter as this term does not affect the shape of the LLH.
Term_Pooledni = sum(log(gamma(n_theta_i+1)),2);

%LLH
LLHnotScaled = exp(Term_Jaz - Term_Pooledfi - Term_Pooledni);

%backup
output.Term_Jaz_meanSim = Term_Jaz;
output.Term_Pooledfi_meanSim = Term_Pooledfi;
output.Term_Pooledni_meanSim = Term_Pooledni;
output.LLHnotScaled = LLHnotScaled;

%plot
if strcmp(Optiondisplay,'display=on')
    
    figure('color','w','position',[224 276 1060 297])
    h(1) = subplot(1,6,1);
    plot(Term_Jaz)
    box off
    axis tight
    title({'Term 1: Jaz.approx of LLH','$\exp{\sum_{i=1}^N (n_i \times \log{f_i(\theta)}))}$'},'interpreter','latex')
    
    h(2) = subplot(1,6,[2 3 4]);
    plot(Term_Pooledfi)
    box off
    axis tight
    title('Term 2: $\exp{(- \sum_{i=1}^N (f_i(\theta)))}$','interpreter','latex')
    xlabel('Motion. direction (degree)')
    
    h(3) = subplot(1,6,5);
    plot(Term_Pooledni)
    box off
    axis tight
    ylim([0.9*min(Term_Pooledni) 1.1*max(Term_Pooledni)])
    title('Term 3: $\exp{(- \sum_{i=1}^N (\log(n_i!)))}$','interpreter','latex')
    
    h(4) = subplot(1,6,6);
    plot(LLHnotScaled)
    box off
    axis tight
    title('LLH (not scaled) term1*term2*term3','interpreter','latex')
    
    SLremoveDeadSpace
end
