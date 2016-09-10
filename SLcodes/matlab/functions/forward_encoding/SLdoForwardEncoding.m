
% SLdoForwardEncoding.m
%
%  author: steeve laquitaine
%    date: 140329
% purpose: run forward encoding modeling on matrix of BOLD (m voxels, n 
%           trials in which different motion directions were presented) for "k" 
%           channels. 
%           We get
%           -a matrix of channels outputs for each of the n directions
%           -the hypothetical channels
%           -the direction space on which they were generated
%           -channels preferred directions
%           -the weights (amplitude) of the k channels for the m voxels
%           
%   %usage: 
%   
%         [Ctrkn,channels,dirSpace,preferredDir,Weightmk] = doForwardEncoding(dataTRAINmn,6,1:1:360,'display=off');
% 
%   description: all motion directions must be repeated to make sure the
%   weights we got from linear regression are correct!

function [Ctrkn,channels,dirSpace,preferredDir,Weightmk] = SLdoForwardEncoding(dataTRAINmn,numchannel,motionDir,OptionDisp)

%draw channels equidistant on the space of motion directions: 
%same-preference-neurons population tuning curve half-wave rectified 
%sinusoid to the power of 5 (become narrower).
dirSpace=de2r(1:1:360,0);
preferredDir=360/numchannel:360/numchannel:360;

%60 degrees channels
channels=nan(360,numchannel);
for i=1:numchannel
    channels(:,i)=drawHalfRectifiedSine(dirSpace,de2r(preferredDir(i),0));
    channels(:,i)=channels(:,i).^5;
end

%plot
if strcmp(OptionDisp,'display=on')==1
    C=linspecer(numchannel);
    figure('color','w')
    set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
    plot(dirSpace,channels,'linesmoothing','on')
    xlabel('Motion directions (degrees)')
    set(gca,'xtick',dirSpace(1:44:360),'xticklabel',ra2d(dirSpace(1:44:360)))
    xlim([0 2*pi])
    box off
end

%make channels' output Ctrkn (k channels, n trials)
Ctrkn=channels(motionDir,:)';

%retrieve channels weights Wmk (m voxels, k channels)
%Linear regression with matrix algebra produces the weights
%W relates the matrix of BOLD responses with the Channels
%Basically Btr(theta)=W.*Ctr(theta)
%Ctrkn'*Ctrkn must be invertible to have a unique solution. Otherwise 
%infinity of solution because of multicollinearity.i.e., columns of 
%Ctrkn are linearly dependent.
%if matrix invertible: 
%rank(Ctrkn*Ctrkn')=size(Ctrkn*Ctrkn')
%det(Ctrkn*Ctrkn') very different from 0
%inv(Ctrkn*Ctrkn')*(Ctrkn*Ctrkn')=identity matrix
%I could not be implement the backslash instead of inv.
if rank(Ctrkn*Ctrkn')~=size(Ctrkn*Ctrkn')
    fprintf('The matrix of output channels is not full rank and thus not invertible. The weights solution are not unique')
else
    fprintf('Good...the matrix of output channels is invertible. Thus the weights are unique solutions')
    Weightmk=dataTRAINmn*(Ctrkn')*inv(Ctrkn*Ctrkn');
end


%draw
figure;
set(gcf,'color','w')
subplot(223)
C=linspecer(size(Weightmk,1));
set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
plot(Weightmk','lineStyle','-','linesmoothing','on')
title('Voxels (colors)','fontsize',14)
set(gca,'xtick',1:numel(preferredDir),'XTickLabel',preferredDir)
xlabel({'Channels by their','preferred direction (degrees)'},'fontsize',14)
ylabel('Channels weights','fontsize',14)

%other way of plotting it
subplot(224)
imagesc(Weightmk)
colormap('hot')
colorbar
xlabel({'Channels by their','preferred direction (degrees)'},'fontsize',14)
ylabel('Voxels','fontsize',14)
set(gca,'xtick',1:numel(preferredDir),'XTickLabel',preferredDir)
title('Weights','fontsize',14)

%draw weighted channels for each voxel
%look at weighted channels for voxel 1 for example.
subplot(2,2,[1 2])
hold all
thisVoxel=1;
C=linspecer(size(channels,2));
set(gca,'NextPlot','replacechildren', 'ColorOrder',C);
plot(1:1:360,Weightmk(thisVoxel(ones(size(channels,1),1)),:).*channels,...
    'linesmoothing','on',...
    'linewidth',2)

xlabel('Motion directions (degrees)','fontsize',14)
set(gca,'fontsize',14)
xlim([0 360])
box off
title('Weighted channels for voxel 1','fontsize',14)
ylabel('Channels outputs (au)','fontsize',14)



