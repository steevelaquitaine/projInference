

%   author: steeve laquitaine
%     date: 140718
%  purpose: from distributions of evidences to likelihood of motion
%          direction 
%reference: Stocker and Simoncelli, 2006, NN and Girshick and Simoncelli.
%

%better way
close 
figure('color','w')
%evidences are the maximum likelihood estimate of neurons population
%responses to stimulus repetitions
clf
subplot(2,2,3)
Pm=vmPdfs(1:1:360,1:1:360,3,'norm');
imagesc(Pm)
xlabel('Stimulus motion direction (deg)')
ylabel('Evidence (deg)')
xlim([0 360])
h=colorbar;
ylabel(h, 'Probability');
title('Evidence densities')

%-------------------------------------------------------
%Evidence distribution over stimulus 180 deg repetitions
%-------------------------------------------------------
subplot(2,2,1)
hold all
PeGivenD180 = Pm(:,180);
plot(PeGivenD180,1:1:360,'k','linesmoothing','on')
plot([0 max(PeGivenD180)],[180 180],'k:')
set(gca,'xticklabel',[])
set(gca,'ytick',[0 180 360],'yticklabel',[0 180 360])
xlim([0 max(PeGivenD180)])
xlabel('Probability')
ylabel('Evidence (deg)')
title('Distribution of evidences for stimulus 180 deg')
box off

%flatten plot
pbaspect([.2 1 1])

%position axis at 180 deg (stimulus direction)
axes = SLgetAxes;
position = get(axes(1),'position');
set(axes(1),'position',[position(1) + 0 position(2:end)])




%----------
%likelihood 
%----------
%likelihood of all stimulus direction given evidence 215
subplot(2,2,4)
hold all
PdirsGivenEvid215 = Pm(215,:);
plot(1:1:360,PdirsGivenEvid215,'k','linesmoothing','on')
plot([215 215],[0 max(PdirsGivenEvid215)],'k:')
box off
set(gca,'yticklabel',[])
set(gca,'xtick',[0 215 360],'xticklabel',[0 215 360])
xlim([0 360])
xlabel('Stimulus motion direction (deg)')
title('Likelihood of directions given evidence 215 deg')

%flatten plot
pbaspect([1 0.2 1])

%position axis at 180 deg (stimulus direction)
axes=SLgetAxes;
position = get(axes(1),'position');
set(axes(1),'position',[position(1)-0.05 position(2)+0.03 position(3:end)])




