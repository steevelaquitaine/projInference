

% slfMRIcheckVoxels.m
%
%
% author: steeve laquitaine
%purpose: visualize voxel goodness of fit and 
%
%
% need to run "slfmriClassify2" and "slfmriStckGetInstancedb.m" before


%% prepro testing: check voxel fits with von mises and determine a threshold r^2
%r-square fit = 0.7 display good qualitative fit. I choose that.
%load voxel params
%hist r2 for each coherence
close all
load d
load VoxelsParams.mat
figure('color','w')
subplot(121); hist(rsquared1,0:0.1:1)
box off; xlim([0 1]); xlabel('R^2'); ylabel('Number of voxels'); title('coh=0.06')
slsetHistColor([.5 .5 .5],'none')
subplot(122); hist(rsquared1,0:0.1:1)
box off; xlim([0 1]); xlabel('R^2'); title('coh=0.12')
    slsetHistColor([.5 .5 .5],'none')

%plot voxel tuning for desired r2 for coh 0.06
%choose a r-square 
testR2 = 0.9;
[~,voxnum] = min(abs(testR2-rsquared1)); rsquared1(voxnum)
figure; st1 = slfmrigetVoxTuning(voxnum,'myRandomDir','myRandomCoh=0.06','mySwitch=2',d,'fminsearch');
vline(225,'k:')
ylim([99 101]); xlim([0 360])
fprintf('# of Instances per direction (x-axis): \n'); disp(st1.count)
fprintf('# std per direction: \n'); disp(st1.std')
fprintf('# mean per direction: \n'); disp(st1.mean')


%% prepro testing: set a high threshold selectivity (vm k param) from coh 0.06
%hist k param for each coherence
%filter voxels that were fit with an r2 of at least 70%
R2thresh = 0.9;
ks1Filt = ks1(rsquared1>=R2thresh);
%filter voxels with high goodness of fit
voxFiltnum = 1:length(ks1); voxFiltnum = voxFiltnum(rsquared1>=R2thresh);
%plot
figure('color','w')
subplot(121); hist(ks1Filt,min(ks1Filt):(max(ks1Filt)-min(ks1Filt))/30:max(ks1Filt))
box off; xlim([min(ks1Filt) max(ks1Filt)]); xlabel('Tuning width (vm k)'); ylabel('# of voxels'); title('coh=0.06')
slsetHistColor([.5 .5 .5],'none')
for i = 1 : length(ks1Filt); ks1Filtdeg(i) = SLKtoStdev(ks1Filt(i)); end
subplot(122); hist(ks1Filtdeg,min(ks1Filtdeg):(max(ks1Filtdeg)-min(ks1Filtdeg))/30:max(ks1Filtdeg))
box off;xlim([min(ks1Filtdeg) max(ks1Filtdeg)]); xlabel('Tuning std (deg)'); ylabel('# of voxels'); title('coh=0.06')
slsetHistColor([.5 .5 .5],'none')

min(ks1)
max(ks1)


%% prepro testing: filter voxels with high goodness of fit and strong selectivity
%select min selectivities in the top 10%
%choose voxels that show at least 1% signal change between less and more preferred
%direction
%width of filtered voxels
kThresh = 0.0049;
ks1FiltR2_K = ks1(rsquared1>=R2thresh & ks1>=kThresh);

%plot each desired tuning
testks1Filt =   0.0070797183655382;
[~,voxnum] = min(abs(testks1Filt - ks1Filt)); 
thisVoxnumFilt = voxFiltnum(voxnum);
figure('color','w'); st1 = slfmrigetVoxTuning(thisVoxnumFilt,'myRandomDir','myRandomCoh=0.06','mySwitch=2',d,'fminsearch');
vline(225,'k:')
ylim([98.5 102]); xlim([0 360])
fprintf('-------------------------\n');
fprintf('tuning width (vm K): \n');disp(ks1(thisVoxnumFilt))
fprintf('tuning width (std,deg): \n');disp(round(SLKtoStdev(ks1(thisVoxnumFilt)))) 
