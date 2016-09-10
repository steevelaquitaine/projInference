

% BayesianBoldMain


%output.rPostk(:,i) = rtotal.*exp(beta*(output.weightLLH(:,i) + output.weightPRIOR))./sum(exp(beta*(output.weightLLH(:,i) + output.weightPRIOR)));



%NO PRIOR
dirTrained = 1:1:360;
neuronIDinVoxLLH = repmat(7,1,1000);
neuronIDinVoxPOST = repmat(250,1,1000);
neuronIDinVoxPrior = repmat(1:15,1,65);

[llhBOLDmn,PostBOLDmn,PriorBOLDmn,tuningPrefinPOST,...
    VoxRespAmpPrior,VoxPrefDirllh,VoxPrefDirPost,neuronIDinVoxLLH,...
    neuronIDinVoxPOST,neuronIDinVoxPrior,output] = simBayesianBold(dirTrained,...
    'priorstd=inf',...
    'priormean=225',...
    'numVoxelsMT=1',...
    'numVoxelsLIP=1',...
    'NeurInaVoxelMT=1000',...
    'NeurInaVoxelLIP=1000',...
    neuronIDinVoxLLH,...
    neuronIDinVoxPOST,...
    neuronIDinVoxPrior,...
    15,...
    'sortVoxels=off',...
    'displayNetwork=off',...
    'displayBOLD=off');

%% sensory neurons firing
%draw
close all
figure(1)
hold all
set(gcf,'color','w','position',[564 674 249 143]);
imagesc(output.sensNeurFIRING)
text(360,7,{'neuron','7'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'Presynaptic','sensory neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[1 15],'k--')
plot([0 360],[7 7],'k--')
set(gca,'fontsize',9)
colorbar
axis tight
axis square
title({'Sensory response (sp)'})

%% sensory weights to readout
figure(2)
hold all
set(gcf,'color','w','position',[811 674 249 143]);
imagesc(output.weightLLH)
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'PostSynaptic','readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
set(gca,'fontsize',9)
colorbar
axis tight
axis square
title({'Sensory presynaptic inputs'})

%% prior weights to readout
fig3 = figure(3);
hold all
set(gcf,'color','w','position',[1048 674 249 143]);
imagesc(output.weightPRIOR(:,ones(1,360)));
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'PostSynaptic','readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
colorbar
set(gca,'fontsize',9)
axis tight
axis square
title({'Prior presynaptic inputs','when NO PRIOR'})

mxWP = max(output.weightPRIOR(:));
minWP = min(output.weightPRIOR(:));

%% readout neurons'responses
%observation: when a prior is applied, the maximum response shifts to
%readout neurons which selectivity is closer to the prior.
figure(4)
hold all
set(gcf,'color','w','position',[1202 673 249 143]);
imagesc(output.rPostk)
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'Readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
plot(1:360,1:360,'k--')
colorbar
axis tight
axis square
title({'Readout neurons firing (sp)','when NO PRIOR'})


%% resulting BOLD
% figure(3)
% hold all
% set(gcf,'color','w','position',[811 674 249 143]);
% plot(llhBOLDmn','k','linesmoothing','on')
% xlabel('Motion directions (deg)')
% box off
% 
% figure(4)
% hold all
% set(gcf,'color','w','position',[1048 674 249 143]);
% plot(PostBOLDmn','k','linesmoothing','on')
% [thmax,argmax] = max(PostBOLDmn);
% plot([argmax argmax],[min(PostBOLDmn) thmax],'--')
% text(argmax,thmax,num2str(argmax))
% xlabel('Motion directions (deg)')
% box off
% 
% figure(5)
% hold all
% set(gcf,'color','w','position',[1202 673 249 143]);
% plot(PriorBOLDmn','k','linesmoothing','on')
% xlabel('Motion directions (deg)')
% box off





%% strong prrior
[llhBOLDmn,PostBOLDmn,PriorBOLDmn,tuningPrefinPOST,...
    VoxRespAmpPrior,VoxPrefDirllh,VoxPrefDirPost,neuronIDinVoxLLH,...
    neuronIDinVoxPOST,neuronIDinVoxPrior,output] = simBayesianBold(dirTrained,...
    'priorstd=0',...
    'priormean=225',...
    'numVoxelsMT=1',...
    'numVoxelsLIP=1',...
    'NeurInaVoxelMT=1000',...
    'NeurInaVoxelLIP=1000',...
    neuronIDinVoxLLH,...
    neuronIDinVoxPOST,...
    neuronIDinVoxPrior,...
    15,...
    'sortVoxels=off',...
    'displayNetwork=off',...
    'displayBOLD=off');

%% sensory neurons firing
%draw
figure(5)
hold all
set(gcf,'color','w','position',[563 482 249 143]);
imagesc(output.sensNeurFIRING)
text(360,7,{'neuron','7'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'Presynaptic','sensory neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[1 15],'k--')
plot([0 360],[7 7],'k--')
set(gca,'fontsize',9)
colorbar
axis tight
axis square
title({'Sensory response (sp)'})

%% sensory weights to readout
figure(6)
hold all
set(gcf,'color','w','position',[811 482 249 143]);
imagesc(output.weightLLH)
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'PostSynaptic','readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
set(gca,'fontsize',9)
colorbar
axis tight
axis square
title({'Sensory presynaptic inputs'})

%% prior weights to readout
figure(7)
hold all
set(gcf,'color','w','position',[1049 481 249 143]);
imagesc(output.weightPRIOR(:,ones(1,360)))
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'PostSynaptic','readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
set(gca,'fontsize',9)
axis tight
axis square
title({'Prior presynaptic inputs','when STRONG PRIOR'})
mxSP = max(output.weightPRIOR(:));
minSP = min(output.weightPRIOR(:));
%range is min and max for strongest prior (std=0)
caxis([-1.7e+3 652])
colorbar

%weak prior graph too same scale
set(fig3,'visible','on')
caxis([min([minWP minSP]) max([mxWP;mxSP])])
colorbar

%% readout neurons'responses
%observation: when a prior is applied, the maximum response shifts to
%readout neurons which selectivity is closer to the prior.
figure(8)
hold all
set(gcf,'color','w','position',[1233 482 249 143]);
imagesc(output.rPostk)
text(360,250,{'neuron','250'},'fontsize',8)
xlabel('Motion direction(deg)','fontsize',11)
ylabel({'Readout neurons','tuning (deg)'},'fontsize',11)
plot([225 225],[0 360],'k--')
plot([0 360],[250 250],'k--')
plot(1:360,1:360,'k--')
colorbar
axis tight
axis square
title({'Readout neurons response (sp)','when STRONG PRIOR'})




% %% resulting BOLD
% figure(3)
% plot(llhBOLDmn','r','linesmoothing','on')
% xlabel('Motion directions (deg)')
% axis tight
% box off
% 
% figure(4)
% hold all
% plot(PostBOLDmn','r','linesmoothing','on')
% [thmax,argmax] = max(PostBOLDmn);
% plot([argmax argmax],[min(PostBOLDmn) thmax],'--')
% text(argmax,thmax,num2str(argmax))
% xlabel('Motion directions (deg)')
% legend('flat prior','pref','strong prior')
% lg=legend('boxoff');
% set(legend,'location','bestoutside')
% axis tight
% box off
% 
% figure(5)
% plot(PriorBOLDmn','r','linesmoothing','on')
% xlabel('Motion directions (deg)')
% axis tight
% box off


