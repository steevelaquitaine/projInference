
function channel=simForwardEncoding
%Create channels
subplot(3,1,1)
title('Channels','fontsize',14,'fontweight','bold')
hold all
numchannel=6;
channel=drawForwardEncoding(numchannel);
ylabel('Relative response','fontsize',14)
%xlabel('Motion directions (degrees)','fontsize',14)


%Weight each channel. Channels weights reflect the number of neurons 
%that contribute to each channel. Then display a motion and measure 
%the six channels activities.
subplot(3,1,2)
title('Channels weighting in example voxel activity','fontsize',14,...
    'fontweight','bold')
hold all
weights=[20 18 25 30 22 15];
weightedchannel=weights(ones(360,1),:).*channel;
plot(weightedchannel,'linesmoothing','on')

%one example motion direction
motiondir=45;
channelResponsetoMotion=weightedchannel(motiondir,:)';
plot([motiondir(ones(numchannel,1)) motiondir(ones(numchannel,1))],...
    channelResponsetoMotion(:,1),'.','color',[.5 .5 .5],'markersize',30)
text(motiondir+5,max(channelResponsetoMotion),{' \leftarrow Response to',...
    'motion direction'},'FontSize',11)

motiondir=225;
channelResponsetoMotion=weightedchannel(motiondir,:)';
plot([motiondir(ones(numchannel,1)) motiondir(ones(numchannel,1))],...
    channelResponsetoMotion(:,1),'.','color',[.5 .5 .5],'markersize',30)
text(motiondir+5,max(channelResponsetoMotion),{' \leftarrow Response to',...
    'motion direction'},'FontSize',11)

ylabel('Weighted response','fontsize',14)
xlim([1 360])
box off


%Convert channels responses to voxel Bold activity
subplot(3,1,3)
voxelBoldtoMotions=sum(weightedchannel,2);
area(voxelBoldtoMotions,'edgecolor','none','facecolor',[0.7 0 0])
text(motiondir,voxelBoldtoMotions(motiondir),...
    {' \downarrow a voxel Bold is the','weighted sum of',...
    'the responses to motion'},'FontSize',11)

xlim([1 360])
ylabel('Voxel Bold','fontsize',14)
xlabel('Motion directions (degrees)','fontsize',14)
xlim([1 360])
box off