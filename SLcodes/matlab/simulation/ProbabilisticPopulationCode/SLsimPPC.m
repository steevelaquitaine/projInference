

%SLsimPPC.m
%
%       $Id: SLsimPPC.m $
%        by: steeve laquitaine
%      date: 140606
%   purpose: Create a Probabilistic Population Code. 
%            Motion directions (e.g., a space of feature) are encoded in 
%            neural tuning responses and likelihood over motion directions 
%            to the readout of an actually presented motion direction.
%           
%     usage:
%
%           output = SLsimPPC(180,20,1,'display=on')
%           output = SLsimPPC(180,100,20,'display=on')
%           output = SLsimPPC(180,2,100,'display=on')
%
%Description:
%
%   - A population of N neurons with tuning curves responses
%     to motion direction
%
%   - Tuning curves are von Mises
%     fi(direction) =  Gain*von mises + base firing;
%
%   - Neurons have preferred directions between 0 and 360º.
%   
%   - Spike count has Poisson variability.
%
%Predictions:    
%
%   - The mode of the likelihood matches the presented direction
%   - When neural gain increases, likelihood becomes stronger (lower std)
%
%
%%references:
%
%    http://homepages.inf.ed.ac.uk/pseries/CCN14/lab4.pdf
%    Jazayeri et al,2006, Nat.Neu
%    Ma et al,2006, Nat.Neu
%    Seung, H.S. & Sompolinsky, H. Simple models for reading neuronal 
%    population codes. Proc. Natl. Acad. Sci. USA 90, 10749?10753(1993).
%
% note: works fine

function output = SLsimPPC(stim,maxResp,nbNeuron,Optiondisplay)

%Poisson noise corrupted tuning responses to all directions
output = SLmakeTuningCurves(nbNeuron,maxResp,1,40);
output = addPoissonNoise(output,Optiondisplay);

%Input motion direction to readout
output.stim = stim;

%Responses to input direction
output = makePopulationResponses(output,Optiondisplay);

%Likelihood of all directions given response
output = readLikelihood(output,Optiondisplay);



%Main functions
%--------------
%Poisson noise 
function output = addPoissonNoise(output,Optiondisplay)

%inputs
f = output.f;
nbNeuron = output.nbNeuron;
diSpace = output.diSpace;

%response range
rSpace = 1 : 1 : 60;

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
pref = output.pref;
PRgivenDi = output.PRgivenDi;

%get population response to a stimulus: sample a spike count response to 
%the stimulus from each neuron
rs = nan(1,output.nbNeuron);
for Ni = 1 : output.nbNeuron
    
    rs(Ni) = randsample(rSpace,1,'true',PRgivenDi(:,stim,Ni));
    
end

%output
output.rs = rs;

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
rs = SLmakeColumn(output.rs);
f = output.f;
diSpace = output.diSpace;

%likelihood (LLH) from Ma & Jazayeri 
%should be same
A = nan(length(diSpace),1);
B = sum(log(gamma(rs+1)));
logL = nan(length(diSpace),1);
LLH_SeungSompo = nan(length(diSpace),1);

for di = 1 : length(diSpace)
    
    %mean response to stimulus
    fdi = f(diSpace(di),:);
    
    %1 - Jazayeri et al.,
    A(di) = sum(fdi);
    logF  = SLmakeColumn(log(fdi));
    logL(di) = rs'*logF - A(di) - B;
    
    %2 - Exact LLH
    %Seung, & Sompolinski, 1993; Jazayeri et al., 2006; Ma et al.2006
    LLH_SeungSompo(di) = prod(exp(-fdi).* fdi.^(rs')./gamma((rs')+1));
end

%proba (should I do that?)
LLH_Jazayeri = exp(logL)/sum(exp(logL));
LLH_SeungSompo = LLH_SeungSompo./sum(LLH_SeungSompo);

%output
output.LLH_SeungSompo  = LLH_SeungSompo;
output.LLH_Jazayeri = LLH_Jazayeri;

%mean LLH
output.meanLLH_Jazayeri = SLmakeMean(diSpace,LLH_Jazayeri);
output.meanLLH_SeungSompo  = SLmakeMean(diSpace,LLH_SeungSompo);

%std LLH
output.stdLLH_Jazayeri = SLmakeStd(diSpace,LLH_Jazayeri);
output.stdLLH_SeungSompo  = SLmakeStd(diSpace,LLH_SeungSompo);

%display
output = display3(Optiondisplay,output);


%nested functions
%--------------
function output = display1(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')==1
    
    %see 5 neurons
    screen = get(0,'ScreenSize');
    figure('position',[286 628 883 180],'color','w')
    
    %examples
    if output.nbNeuron >= 5
        ex = output.nbNeuron/5 : output.nbNeuron/5 : output.nbNeuron;
    else
        ex = 1;
    end
    
    output.nex = length(ex);
    for i = 1 : output.nex
        
        %axis
        subplot(1,output.nex+2,i)
        
        %plot
        imagesc(output.PRgivenDi(:,:,ex(i)))
        
        %graphics
        set(gca,'YDir','normal');
        
        %label 1st plot
        if i==1
            ylabel('Spike counts','fontsize',12)
        end
        
        %legend last plot
        if i==round(output.nex/2)
            SzD1 = size(output.PRgivenDi,2);
            SzD2 = size(output.PRgivenDi,1);
            text(SzD1/6,SzD2/1,'P(responses)')
            xlabel('Mot. dir.(º)')
        end
        axis tight
        box off
    end
    
    %color
    map = colormap('gray');
    colormap(flipud(map));
end

%population response sorted by neurons preference
function output = display2(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')==1
    
    subplot(1,output.nex+2,output.nex+1)
    axis square
    
    %responses
    scatter(output.pref,output.rs,40,'markerfacecolor','b')
    
    %labels
    set(gca,'xtick',0:90:360,'xticklabel',0:90:360)
    xlim([1 360])
    ylim([0 max(max(output.rs))+1.1])
    
    xlabel('Neur. (by pref. (º))')
    SzD1 = size(output.pref,2);
    SzD2 = size(output.rs,1);
    text(SzD1/2,max(max(output.rs))+1.1,'Pop.resp.')
    
end

%likelihood
function output = display3(Optiondisplay,output)

if strcmp(Optiondisplay,'display=on')
    
    %draw
    subplot(1,7,7)
    axis square
    hold all
    
    %plot(logL)
    plot(output.LLH_SeungSompo,'-g','linesmoothing','on','linewidth',2)
    plot(output.LLH_Jazayeri,'--','color',[0.9 0.2 0.2],'linesmoothing','on','linewidth',2)
    
    %labels
    legend('Exact LLH (Seung,Sompo)','Jazayeri LLH')
    legend('boxoff')
    
    set(gca,'xtick',0:90:360,'xticklabel',0:90:360)
    xlabel('Mot. dir (º)')
    ylabel('Likelihood')
    ylim([0 max(max([output.LLH_SeungSompo; output.LLH_Jazayeri]))])
    xlim([0 360])
    box off
    SLremoveDeadSpace

end