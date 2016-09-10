

%SLsimPPCafterAdaptation.m
%
%       $Id: SLsimPPCafterAdaptation.m $
%        by: steeve laquitaine
%      date: 141107
%
%   purpose: Create a Probabilistic Population Code. 
%            Motion directions (e.g., a space of feature) are encoded in 
%            neural tuning responses and likelihood over motion directions 
%            to the readout of an actually presented motion direction.
%            Tuning responses are adapted.
%           
%     usage:
%
%           output = SLsimPPCafterAdaptation(180,1,100,20,225,0,0,[0:1:60],'display=on')
%        
% 
%Description:
%
%   - A population of N neurons with von mises tuning to direction.
%
%   - Tuning curves are von Mises
%     fi(direction) =  Gain*von mises + base firing;
%
%   - Neurons have preferred directions between 0 and 360º.
%   
%   - Spike count has Poisson variability.
%
%   - Tuning curves are adapted at 225º(adapter) by a circular multiplicative 
%     gain of 28.3º std (k=4.7) with amplitude 0.4 (see Benucci et al.,).
%
%
%Predictions:    
%
%     Strong motion (maxResp = 60 spikes)
%       output = SLsimPPCafterAdaptation(180,1,60,100,225,4.7,1,[1:1:61],'display=on')
%       - LLH peaks at input direction
%       - LLH is strong (low std)
%
%     Weak motion (maxResp = 2 spikes)
%       output = SLsimPPCafterAdaptation(180,1,2,100,225,4.7,1,[1:1:61],'display=on')
%       - LLH shows two peaks: at input direction and at adapter
%
%      WHY ?
%      
%References:
%
%    PPC
%       http://homepages.inf.ed.ac.uk/pseries/CCN14/lab4.pdf
%       Jazayeri et al,2006, Nat.Neu
%       Ma et al,2006, Nat.Neu
%       Seung, H.S. & Sompolinsky, H. Simple models for reading neuronal 
%       population codes. Proc. Natl. Acad. Sci. USA 90, 10749?10753(1993).
%
%    Adaptation
%       Benucci, A., Saleem, A. B. & Carandini, M. 
%       Adaptation maintains population homeostasis in primary 
%       visual cortex. Nature Publishing Group 16, 724?729 (2013).
%
% note: works fine

function output = SLsimPPCafterAdaptation(stim,base,maxResp,nbNeuron,Adap,kAdap,magAdap,...
    rSpace,Optiondisplay)

%Tuning responses to direction space
output = SLmakeTuningCurves(nbNeuron,maxResp,base,40);


%ParametePopResp
output.stim = stim;
output.Adap = Adap;
output.kAdap = kAdap;
output.magAdap = magAdap;
output.rSpace = rSpace;

%Adapt to a repeated stimulus (Adapter,e.g.,225º)
output = AdaptTuningCurves(output,Optiondisplay);

%Poisson noise
output = addPoissonNoise(output,Optiondisplay);

%Responses to input direction
output = makePopulationResponses(output,Optiondisplay);

%Likelihood of all directions given response
output = readLikelihood(output,Optiondisplay);



%Main functions
%--------------
%Adaptation of tuning curves
function output = AdaptTuningCurves(output,Optiondisplay)

%inputs
Adap    = output.Adap;
kAdap   = output.kAdap;
magAdap = output.magAdap;
nbNeuron= output.nbNeuron;
diSpace = output.diSpace;
D       = output.D;
f       = output.f;
pref    = output.pref;

%adaptation gain
AdapGain = 1 - magAdap.*vmPdfs(diSpace,Adap,kAdap,'NotNorm');

%warning
if any(AdapGain < 0)
    sprintf(['(SLadaptTuningCurves) Adaptation gain is < 0. You should',...
        ' change adaptation parameters to make the gain > 0.'])
    keyboard
end

posAdaptedPref = nan(nbNeuron,1);
for Ni = 1 : nbNeuron
    posAdaptedPref(Ni) = find(diSpace==pref(Ni));
end
AdapGainByPref = AdapGain(posAdaptedPref)';

%Non adapted tuning curves becomes adapted
output.f_nonAdapted = f;
output.f = AdapGainByPref(ones(D,1),:).*f;

%output
output.AdapGain = AdapGain;
output.AdapGainByPref = AdapGainByPref;
output.AdaptedPref = diSpace(posAdaptedPref');

%display
display0(Optiondisplay,output)

%Poisson noise 
function output = addPoissonNoise(output,Optiondisplay)

%inputs
f = output.f;
nbNeuron = output.nbNeuron;
diSpace = output.diSpace;

% response range (spike count)
% must be defined as a Natural number N = {0:1:+inf}
rSpace = output.rSpace;

%P(responses|direction,neuron)
%Probability are defined by the Poisson noise
%note: for non integer values the gamma function interpolates
%factorial function. Factorial(n) = gamma(n+1);
PRgivenDi = nan(numel(rSpace),length(diSpace),nbNeuron);
for Di = 1 : length(diSpace)
    for Ni = 1 : nbNeuron
        
        PRgivenDi(:,Di,Ni) = (exp(-f(Di,Ni)) * f(Di,Ni).^rSpace)./gamma(rSpace+1);
        
        %sum to 1
        PRgivenDi(:,Di,Ni) = PRgivenDi(:,Di,Ni)/sum(PRgivenDi(:,Di,Ni));
    end
end

%output
output.nbNeuron = nbNeuron;
output.PRgivenDi = PRgivenDi;
output.rSpace = rSpace;

%display
output = display1(Optiondisplay,output);

%Population responses
function output = makePopulationResponses(output,Optiondisplay)

%inputs
stim = output.stim;
rSpace = output.rSpace;
PRgivenDi = output.PRgivenDi;

%get population response to a stimulus: sample a spike count response to 
%the stimulus from each neuron
PopResp = nan(1,output.nbNeuron);
for Ni = 1 : output.nbNeuron
    
    PopResp(Ni) = randsample(rSpace,1,'true',PRgivenDi(:,stim,Ni));
    
end

%output
output.PopResp = PopResp;

%diplay
output = display2(Optiondisplay,output);

%read 
function output = readLikelihood(output,Optiondisplay)
%logLLH(stimuli|response) is 
%   Sum(neurons'responses to stim weighted by log(mean response for this
%   stim))
%
%where the mean response for the stim is defined by the tuning curve
%
%note: can we really safely discard the 2nd
%term in the equation as stated by Jazayeri et al.?
%The second term can be discarded when the 
%representation is homogeneous. Is it the case?

%inputs
PopResp = SLmakeColumn(output.PopResp);
f = output.f;
diSpace = output.diSpace;

%likelihood (LLH) from Ma & Jazayeri 
%should be same
A = nan(length(diSpace),1);
B = sum(log(gamma(PopResp+1)));
logL = nan(length(diSpace),1);
%LLH_SeungSompo = nan(length(diSpace),1);

for di = 1 : length(diSpace)
    
    %tuning response
    fdi = f(diSpace(di),:);
    
    %1 - log(LLH) Jazayeri et al.,2006 formulation
    A(di) = sum(fdi);
    logF  = SLmakeColumn(log(fdi));
    logL(di) = PopResp'*logF - A(di) - B;

    %2 - Exact LLH
    %Seung, & Sompolinski, 1993; Jazayeri et al., 2006; Ma et al.2006
    %numerical limits for large base of max responses change of the neurons
    %or for large number of neurons
    %LLH_SeungSompo(di) = prod(exp(-fdi).* fdi.^(PopResp')./gamma((PopResp')+1));
end

%Proba (should I do that?)
%numerical limits for large base of max responses change of the neurons
%of for large number of neurons
LLH_Jazayeri = exp(logL)/sum(exp(logL));
%LLH_SeungSompo = LLH_SeungSompo./sum(LLH_SeungSompo);

%output
output.logLLH_Jazayeri = logL;
%output.LLH_SeungSompo  = LLH_SeungSompo;
output.LLH_Jazayeri = LLH_Jazayeri;

%mean LLH
[~,posModeJaz] = max(LLH_Jazayeri);
%[~,posModeMa]  = max(LLH_SeungSompo);
output.modeLLH_Jazayeri = diSpace(posModeJaz);
%output.modeLLH_SeungSompo = diSpace(posModeMa);

%std LLH
output.stdLLH_Jazayeri = SLmakeStd(diSpace,LLH_Jazayeri);
%output.stdLLH_SeungSompo  = SLmakeStd(diSpace,LLH_SeungSompo);

%display
output = display3(Optiondisplay,output);





%nested functions
%----------------
%Adaptation gain and Tuning curves
function output = display0(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')==1
    
    %figure
    figure('position',[151 676 411 387],'color','w')
        
    %plot
    %Adaptation gain
    %---------------
    subplot(3,1,1)
    plot(output.AdaptedPref,output.AdapGainByPref,'.k-',...
        'Markerfacecolor','k')
    box off
    
    %labels
    set(gca,'xtick',output.AdaptedPref(1:36:end),'xticklabel',...
        output.AdaptedPref(1:36:end))
    xlim([1 360])
    xlabel('Neurons by pref (º)')
    ylabel('Adaptation (X gain)')
    
    %10 example adapted tuning curves
    %--------------------------------
    h = subplot(3,1,2);
    
    %color neurons
    C = SLlinspecer(10);
    set(h,'NextPlot','replacechildren','ColorOrder',C)
    
    %plot
    if output.nbNeuron < 10
        plot(output.diSpace,output.f,'linesmoothing','on')
        box off
    else 
        expos = round(1 : output.nbNeuron/10 : output.nbNeuron);
        plot(output.diSpace,output.f(:,expos),'linesmoothing','on')
    end
    %set(gca,'ColorOrder',C);

    
    %labels
    set(gca,'xtick',output.AdaptedPref(expos),'xticklabel',...
        output.AdaptedPref(expos))
    xlim([1 360])
    ylim([0 max(output.f(:))])
    xlabel('Mot.dir.(º)')
    ylabel('Resp.(spikes)')

    %All adapted tuning curves
    %-------------------------
    subplot(3,1,3)
    mesh(output.f)
    view(gca,[-24.5 24]);
%     SLremoveDeadSpace
    colorbar
end

%Example neuron with added Poisson noise
function output = display1(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')==1
    
    %example neurons
    figure('position', [155 425 883 180],'color','w')
    ex = round(output.nbNeuron/2);
    
    output.nex = length(ex);
    for i = 1 : output.nex
        
        %axis
        subplot(1,output.nex+2,i)
        
        %plot
        imagesc(output.PRgivenDi(:,:,ex(i)))
        
        %this plots the actual average response
        %to each motion direction (tuning curves) after 
        %addition of Poisson noise
        %plot(sum(output.PRgivenDi(:,:,ex(i)).*repmat(rSpace',1,output.nbNeuron)))
        
        %graphics
        set(gca,'YDir','normal');
        
        %label 1st plot
        if i==1
            ylabel('Spike counts','fontsize',12)
        end
        
        %legend last plot
        SzD1 = size(output.PRgivenDi,2);
        SzD2 = 3*output.maxResp;
        text(SzD1/6,SzD2/1,'A noisy neuron')
        xlabel('Mot. dir.(º)')
        axis tight
        box off
        ylim([0 3*output.maxResp])
    end
    set(gca,'xtick',0:90:360,'xticklabel',0:90:360)
    
    %color
    map = colormap('gray');
    colormap(flipud(map));
end

%Population response sorted by neurons preference
function output = display2(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')==1
    
    subplot(1,output.nex+2,output.nex+1)
    axis square
    
    %responses
    scatter(output.pref,output.PopResp,20,'markerfacecolor','b')
    
    %labels
    set(gca,'xtick',output.pref(1:36:end),'xticklabel',output.pref(1:36:end))
    xlim([1 360])
    ylim([0 max(max(output.PopResp))+1.1])
    
    xlabel('Neur. (by pref. (º))')
    SzD1 = size(output.pref,2);
    SzD2 = size(output.PopResp,1);
    text(SzD1/2,max(max(output.PopResp))+1.1,'Pop.resp.')
    
end

%likelihood
function output = display3(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')
    
    %draw
    subplot(1,7,7)
    hold all
    maxLLH = max(output.LLH_Jazayeri);
    
    %LLH
    %plot(output.diSpace,output.LLH_SeungSompo,'-g','linewidth',1.5,...
    %    'linesmoothing','on')
    plot(output.diSpace,output.LLH_Jazayeri,'--','color',[0.9 0.2 0.2],...
        'linewidth',1.5,'linesmoothing','on')
    
    text(output.modeLLH_Jazayeri-90,maxLLH,[num2str(round(output.modeLLH_Jazayeri)),' deg'],...
        'color','r','fontsize',8)
    
    %readout
    plot([output.modeLLH_Jazayeri output.modeLLH_Jazayeri],[0 maxLLH],...
        'r:','linewidth',1)
    
    %adapter
    plot([output.Adap output.Adap],[0 maxLLH],'b:','linewidth',1,...
        'linesmoothing','on')
    
    text(output.Adap,maxLLH,['Adapter is ',num2str(output.Adap),' deg'],...
        'color','b','fontsize',8)
    
    %labels
    %legend('Exact LLH (Seung,Sompo)','Jazayeri LLH')
    legend('Jazayeri LLH')
    legend('boxoff')
    set(gca,'xtick',0:90:360,'xticklabel',0:90:360)
    xlabel('Mot. dir (º)')
    ylabel('Likelihood')
    %ylim([0 max(max([output.LLH_SeungSompo; output.LLH_Jazayeri]))])
    ylim([0 max(max(output.LLH_Jazayeri))])
    xlim([0 360])
    box off
    SLremoveDeadSpace
end