
%convHRF.m
%
% author: steeve Laquitaine
%   date: 140422 last modification 140426
%purpose: convolve a time series (stimulus presentation) with the 
%         hemodynamic response function.
%
%  usage: 
%
%          
%          tsteps         = 0.1;
%          stimDuration   = 3.2;
%          trialDuration  = 16;
%          numTrials      = 1;
%          tSeries        = zeros(1,round(trialDuration./tsteps))';
%          tSeries(1 : round(stimDuration./tsteps)) = 1;
%          tSeries        = repmat(tSeries,numTrials,1);
%          [tSeriesConvWithHRF,hrf] = convHRF(tSeries,tsteps,'OptionDisp=on');
%
%
%
% Description: Convolution is :
%              w(k) = sum_j( u(j) * v(k-j+1) )
%              or (f*g)(t) = sum_0_to_t(f(tau)*g(t-tau)) d_tau
%              
%              time reverse g(t) (i.e., t --> -t)
%              shift by t
%              then sum intersecting area between f and g.
%              
% reference: https://labs.psych.ucsb.edu/ashby/gregory/Matlab.html
% 
% notes: 
% - convolved BOLD increases and saturates when stimulus lasts.
%   because of the temporal properties of the HRF. Summing BOLD in early phase
%   of HRF (before it peaks) changes the ouput BOLD a lot particularly  
%   closer to HRF peak that adds the most. 
%
% - Then later phase of the HRF go down
%   slowly which causes output BOLD to saturate (steady phase) because. If
%   it dropped suddenly, summed output BOLD would go down too. It's like a
%   leaky balloon. Need to think more about it.


function [tSeriesConvWithHRF,hrf,tsteps] = convHRF(tSeries,tsteps,OptionDisp)
 
%time axis
t = 0 : tsteps : 30;
 
%number of time series to convolve
numtSeries = size(tSeries,2);

%number of time steps
lengthSeries = size(tSeries,1);

%hrf
T0 = 0;n = 4;lamda = 2;
hrf = ((t-T0).^(n-1)).*exp(-(t-T0)/lamda)/((lamda^n)*factorial(n-1));
hrf = hrf/sum(hrf);

%%time series convolved with HRF
for j = 1 : numtSeries
    tSeriesConvWithHRF(:,j)=conv(tSeries(:,j),hrf);
end

%keep the signal period corresponding to the stimulus time series period
tSeriesConvWithHRF=tSeriesConvWithHRF(1:lengthSeries,:);

%plot the first time series as an example
if strcmp(OptionDisp,'OptionDisp=on')
    
    %stimulus
    figure('color','w')
    subplot(2,1,1)
    bar(tSeries(:,1),'k')
    xlim([0 lengthSeries])
    set(gca,'xtick',1:10:lengthSeries,'xticklabel',t*10,'fontsize',14)
    box off
    title('Stimulus time series')

    
    %BOLD
    subplot(2,1,2)
    plot(tSeriesConvWithHRF(:,1),'linewidth',2,'linesmoothing','on','color',[.7 0 0])
    box off
    xlim([0 lengthSeries])
    set(gca,'xtick',1:10:lengthSeries,'xticklabel',t*10,'fontsize',14)
    ylabel('BOLD signal (au)','fontsize',14)
    xlabel('Time (seconds)','fontsize',14)
    box off
    title('Convolution of time series and hrf')

end