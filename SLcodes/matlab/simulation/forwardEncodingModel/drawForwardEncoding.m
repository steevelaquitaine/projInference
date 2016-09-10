
%author: Steeve Laquitaine
  %date: 140117
  %purpose: draw forward encoding predictions of activity of an fMRI voxel.
    %usage: channel=drawForwardEncoding(6)  
    %to create six equidistant channels
  %reference: Kok, P., Brouwer, G. J., van Gerven, M. A. J. & de Lange, F. P. Prior expectations bias sensory representations in visual cortex. Journal of Neuroscience 33, 16275?16284 (2013).

function [channel,preferredDir]=drawForwardEncoding(numchannel,OptionDisp)
%draw channels equidistant on the space of motion directions: 
%same-preference-neurons population tuning curve half-wave rectified 
%sinusoid to the power of 5 (become narrower).
motiondir=de2r(1:1:360,0);
preferredDir=360/numchannel:360/numchannel:360;

%60 degrees selective channel
channel=nan(360,numchannel);
for i=1:numchannel
    channel(:,i)=drawHalfRectifiedSine(motiondir,de2r(preferredDir(i),0));
    channel(:,i)=channel(:,i).^5;
end

if strcmp(OptionDisp,'display=on')==1
    C=linspecer(numchannel);
    set(gcf,'color','w')
    set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
    plot(motiondir,channel,'linesmoothing','on')
    xlabel('Motion directions (degrees)')
    %graphics
    set(gca,'xtick',motiondir(1:44:360),'xticklabel',ra2d(motiondir(1:44:360)))
    xlim([0 2*pi])
    box off
end