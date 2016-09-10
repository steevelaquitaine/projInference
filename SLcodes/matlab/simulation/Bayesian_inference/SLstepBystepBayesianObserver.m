

%SLStepByStepIdealObserver


%TRIAL 1
%% Frame 1
%Sensory response is not noisy. So a sensory response always indicates one
%and only one motion direction. 
clf
set(gcf,'color','w')
hold all

Motdir = 1:1:360;
llh = zeros(360,1);
llh(135) = 1;
bar(Motdir,llh,'facecolor','k')
set(gca,'xtick',10:25:355,'xticklabel',10:25:355,'fontsize',8)
ylabel('Probability')
xlabel('Motion direction (º)')
box off
SLConventionUp(gcf)

%% Frame 2
%Sensory response is noisy. So a sensory response can many different motion
%directions are likely given this sensory response.
%TRIAL 1
%figure;
clf
hold all

%plot llh
llh = vmPdfs(Motdir,135,2,'norm');
area(Motdir,llh,'facecolor','k','edge','none')
plot([135 135],[0 max(llh)],'w')
prior = vmPdfs(Motdir,225,2.2,'norm');
posterior = llh.*prior/sum(llh.*prior);
ylim([0 max([llh; prior;posterior])])


%graphics
set(gca,'xtick',10:25:355,'xticklabel',10:25:355,'fontsize',8)
box off
SLConventionUp(gcf)


%% Frame 3
%plot prior
area(Motdir,prior,'facecolor',[.8 0 0],'edge','none')
plot([225 225],[0 max(prior)],'w')

%% Frame 4 
%plot posterior
area(Motdir,posterior,'facecolor',[.9 .5 .5],'edge','none')
meanpos = sum(posterior.*Motdir');
plot([meanpos meanpos],[0 max(posterior)],'w')

%% Frame 5
%TRIAL 2
%The same motion direction 135º is repeated in a different trial, trial 2
%and produces sensory response 2 that is most often produced by another
%motion direction, 180º for example.
MostLikelyDirOfSensoryResp2 = 90;

%plot llh trial 2
llh2 = vmPdfs(Motdir,MostLikelyDirOfSensoryResp2,2,'norm');
area(Motdir,llh2,'facecolor',[.5 .5 .5],'edge','none')
plot([MostLikelyDirOfSensoryResp2 MostLikelyDirOfSensoryResp2],[0 max(llh2)],'w')

%% Frame 6
posterior2 = llh2.*prior/sum(llh2.*prior);
area(Motdir,posterior2,'facecolor',[1 .7 .7],'edge','none')
meanpos2 = sum(posterior2.*Motdir');
plot([meanpos2 meanpos2],[0 max(posterior2)],'w')






%TRIAL 2