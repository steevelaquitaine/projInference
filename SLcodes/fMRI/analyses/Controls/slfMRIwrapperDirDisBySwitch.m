

%slfMRIwrapperDirDisBySwitch.m
%
%
% author: steeve laquitaine
%   date: 160310
%purpose: plot distribution of stimulus motion directions sorted by 
%         switching variable
%
%Description: 
%   1) initialize data and analysis info
%   2) make database
%   3) plot distribution of directions sorted by switching variable

%initialize analysis and make fMRI database
slfmriInitAnalysisTaskDotDirfMRI05
[d,o] = slfMRImakeDatabase(o,varargin);

%% plot distribution of directions sorted by switching variable
figure('color','w'); 
subplot(1,2,1)
hold all
h1 = hist(d.myRandomDir(d.mySwitch==1 & d.myRandomCoh==0.06),1:1:360);
bar([1:1:360]-4,h1,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h2 = hist(d.myRandomDir(d.mySwitch==2 & d.myRandomCoh==0.06),1:1:360);
bar([1:1:360]+4,h2,10,'facecolor','k','lineStyle','none')
ylabel('Count')
xlabel('Motion direction (deg)')
xlim([0 360])
legend('Sw-p','Sw-d')
alpha(0.9)
set(gca,'xtick',unique(d.myRandomDir),'xticklabel',unique(d.myRandomDir))
vline(225,':k')
title('Motion direction dist. by switching for 6% coh')
box off

subplot(1,2,2)
hold all
h3 = hist(d.myRandomDir(d.mySwitch==1 & d.myRandomCoh==0.12),1:1:360);
bar([1:1:360]-4,h3,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h4 = hist(d.myRandomDir(d.mySwitch==2 & d.myRandomCoh==0.12),1:1:360);
bar([1:1:360]+4,h4,10,'facecolor','k','lineStyle','none')
ylabel('Count')
xlabel('Motion direction (deg)')
xlim([0 360])
legend('Sw-p','Sw-d')
alpha(0.9)
set(gca,'xtick',unique(d.myRandomDir),'xticklabel',unique(d.myRandomDir))
vline(225,':k')
title('For 12% coh')
box off
