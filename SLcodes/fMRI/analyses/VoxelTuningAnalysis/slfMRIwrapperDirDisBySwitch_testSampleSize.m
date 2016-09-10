

%slfMRIwrapperDirDisBySwitch_testSampleSize.m
%
%
% author: steeve laquitaine
%   date: 160423
%purpose: plot distribution of directions sorted by switching variable
%
%Description: 
%   1) initialize data and analysis info
%   2) make database (take half the samples)
%   3) plot distribution of directions sorted by switching variable

%initialize analysis and make fMRI database
slfmriInitAnalysisTaskDotDirfMRI05
[d,o] = slfMRImakeDatabase(o,varargin);

%% cut # of trials by half
%backup d
dcut = d;

%find trials each direction
ix15_06 = find(d.myRandomDir==15 &d.myRandomCoh==0.06);
ix85_06 = find(d.myRandomDir==85 &d.myRandomCoh==0.06);
ix155_06 = find(d.myRandomDir==155 &d.myRandomCoh==0.06);
ix225_06 = find(d.myRandomDir==225 &d.myRandomCoh==0.06);
ix295_06 = find(d.myRandomDir==295 &d.myRandomCoh==0.06);

%reduce to half the trials
ix15_06half = ix15_06(1:floor(length(ix15_06)/2));
ix85_06half = ix85_06(1:floor(length(ix85_06)/2));
ix155_06half = ix155_06(1:floor(length(ix155_06)/2));
ix225_06half = ix225_06(1:floor(length(ix225_06)/2));
ix295_06half = ix295_06(1:floor(length(ix295_06)/2));
ixhalf = [ix15_06half;ix85_06half;ix155_06half;ix225_06half;ix295_06half];

%update the database
dcut.instances = dcut.instances(ixhalf,:);
dcut.stimvols = dcut.stimvols(ixhalf,:);
dcut.myRandomCoh = dcut.myRandomCoh(ixhalf,:);
dcut.myRandomDir = dcut.myRandomDir(ixhalf,:);
dcut.mySwitch = dcut.mySwitch(ixhalf,:);


%% plot distribution of directions sorted by switching variable
figure('color','w'); 
subplot(1,2,1)
hold all
h1 = hist(dcut.myRandomDir(dcut.mySwitch==1 & dcut.myRandomCoh==0.06),1:1:360);
bar([1:1:360]-4,h1,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h2 = hist(dcut.myRandomDir(dcut.mySwitch==2 & dcut.myRandomCoh==0.06),1:1:360);
bar([1:1:360]+4,h2,10,'facecolor','k','lineStyle','none')
ylabel('Count')
xlabel('Motion direction (deg)')
xlim([0 360])
legend('Sw-p','Sw-d')
alpha(0.9)
set(gca,'xtick',unique(dcut.myRandomDir),'xticklabel',unique(dcut.myRandomDir))
vline(225,':k')
title('Motion direction dist. by switching for 6% coh')
box off

subplot(1,2,2)
hold all
h3 = hist(dcut.myRandomDir(dcut.mySwitch==1 & dcut.myRandomCoh==0.12),1:1:360);
bar([1:1:360]-4,h3,10,'facecolor',[0 0.5 .7],'lineStyle','none')
h4 = hist(dcut.myRandomDir(dcut.mySwitch==2 & dcut.myRandomCoh==0.12),1:1:360);
bar([1:1:360]+4,h4,10,'facecolor','k','lineStyle','none')
ylabel('Count')
xlabel('Motion direction (deg)')
xlim([0 360])
legend('Sw-p','Sw-d')
alpha(0.9)
set(gca,'xtick',unique(dcut.myRandomDir),'xticklabel',unique(dcut.myRandomDir))
vline(225,':k')
title('For 12% coh')
box off
