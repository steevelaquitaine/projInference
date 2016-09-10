
     %Author: Steeve Laquitaine
       %date: 140428 - last modification: 140510
  %file name: SLfMRIeventRelatedAnalysis
  
      %usage: The code runs step by step. 
      
%purpose: Run event-related analysis on the main task experiment. This
%areas whose activity is modulated by motion coherence.
 
%list of analyses:
%- display a map of r2 values calculated from fitting hdrs with time Series 
%of BOLD percent signal changes on a flat map.
%- bar number of voxels with r2 greater than increasing r2 thresholds
%- draw event-related plot for each voxel sorted by their r2.


%backup command window
diary Log_SLfMRIeventRelatedAnalysis_140505.txt

%backup m.file with date and time
filename=matlab.desktop.editor.getActiveFilename;
SLBackup(filename)

%Run correlaton analysis on the localizer experiment.
cd /Users/steeve/fMRI_data/s02520140501/
%cd /Users/steeve/fMRI_data/s02520140425/
%cd /Users/steeve/fMRI_data/s02520140220/

%launch view
mrLoadRet([])

%set the task that we want to analyse
taskNum=1;
%taskNum=2;
phaseNum=1;
segmentNum=2;

%set scans of the task we want to analyse
scanNumFromDescription='main task';

%check if a concatenation group exist
groupNames=viewGet(getMLRView,'groupNames');

%if it does not exist, concatenate
if sum(strcmp(groupNames,'Concatenation'))~=1;

    %set motion compensation group for concatenation
    curGroup='MotionComp';
    viewSet(getMLRView,'curGroup',curGroup);
    scanNum=viewGet(getMLRView,'scanNumFromDescription',...
        scanNumFromDescription,curGroup);
    
    %detrend and high-pass filter the scans you select, get rid of any
    %junk frames, do a de-trending step using a hi-pass filter if you
    %specify that, and then concatenate the resulting time series into a
    %single time series that will be saved into a new group called
    %Concatenation. The mat files are saved in folder "concatenation".
    %set concatenation parameters
    [~,params]=concatTSeries(getMLRView,[],'justGetParams=1',...
        'defaultParams=1',...
        'scanList',scanNum);
    
    %concatenate (stack preprocessed scans into one long scan)
    concatTSeries(getMLRView,params);
    fprintf('\n %12s \n','(SLfMRIeventRelatedAnalysis) Motion compensation scans have been concatenated')    
end

%set concatenated group and scan for event-related analysis
viewSet(getMLRView,'curGroup','Concatenation');
nScans=viewGet(getMLRView,'nScans');
viewSet(getMLRView,'curScan',nScans);

%load flat map anatomy (requires 2 files e.g., 's025_leftFlat.hdr/img/mat')
%some flat maps have been created for retinotopy. We can just use those
%here too. We use those made for the first pilot.
%left
loadAnat(getMLRView,'s025_leftflat.hdr',...
    '/Users/steeve/fMRI_data/s02520140220/Anatomy');

%%Load and overlay ROIs
%myROIname='V1toMT';
myROIname='allMT';
     
%left
loadROI(getMLRView,{'hV4','V3v','V2v','V1','V2d','V3d',...
        'V3a','V3b','LO1','LO2','allMT',myROIname},0,...
              '/Users/steeve/fMRI_data/s02520140220/ROIs')
          
%show ROIs
viewSet(getMLRView,'showrois','all perimeter');

%color ROIs perimeter
roiNum=viewGet(getMLRView,'roiGroup');
for j=1:numel(roiNum)
    viewSet(getMLRView,'roiColor','white',roiNum(j));
end


%Event-related r2 map
%---------------------
%get the acquired volumes in our task (taskNum=1, phaseNum=1) for our 
%variable at motion stimulus onset (segmentNum=2)
%run event-related analysis
%set parameters
[~,params]=eventRelated(getMLRView,[],'justGetParams=1',...
    'defaultParams=1','scanList',nScans);

%apply filtering
params.applyFiltering=1;

%set my variable for the (unique or most recent) concatenation scan
myVariable='myRandomCoh';
params.scanParams{nScans}.varname=myVariable;

%set the task
params.scanParams{nScans}.taskNum=taskNum;
params.scanParams{nScans}.segmentNum=segmentNum;

%display an r2 map from the event-related analysis
v=eventRelated(getMLRView,params);


%event-related plot
%------------------
%make an event-related plot based on ROI r2's thresholded map. We keep the
%time series of the voxels that show the greater r2 values (voxels activity 
%that are modeled best). We average those timeSeries into on time series
%and we use this time series to estimate an hemodynamic response function
%for our different stimuli (e.g., coherences).

%make an average time Series for all the voxels that have a r2 greater than
%a certain threshold
%load ROI
myROI=loadROITSeries(getMLRView,myROIname);

%converting the scan coordinates of the ROI to linear coordinates
scanDims=viewGet(getMLRView,'scanDims');
myROI.linearScanCoords=sub2ind(scanDims,myROI.scanCoords(1,:),...
    myROI.scanCoords(2,:),myROI.scanCoords(3,:));

%get r2 map that will permit to choose greatest r2 voxels
r2=viewGet(getMLRView,'overlayData',nScans);

%select ROI's r2s
myROI.r2=r2(myROI.linearScanCoords);

%average time series of the greatest r2 voxels
myr2Thresh=0.2;
tSeries=mean(myROI.tSeries(myROI.r2>myr2Thresh,:));

%number of voxels with r2 greater than the threshold
numVoxelsTotal=sum(myROI.r2>0);
numVoxels=sum(myROI.r2>myr2Thresh);

%choose the number of hemodynamic response function to deconvolve (i.e., 
%number of stimuli) from the average time Series.
stimvol=getStimvol(getMLRView,myVariable,'taskNum',taskNum,'phaseNum',...
    phaseNum,'segmentNum',segmentNum);
nhdr=length(stimvol);

%choose the length of the hemodynamic response functions in blocks
hdrlen=5;

%make a stimulus convolution matrix (instances in rows, hdr's time X number 
%of stimuli in column. Values are 1 at the times when we estimate hdr
%values).
scm=makescm(getMLRView,hdrlen,1,stimvol);

%deconvolve the hdrs from the mean time series
d=getr2timecourse(tSeries,nhdr,hdrlen,scm,viewGet(getMLRView,'framePeriod'));

%plot the result (event-related plots).
figure('color','w');
hold all
ercolors=[0 0 0; 1 0 0];
for j=1:nhdr
%     errorbar(d.time,d.ehdr(j,:),d.ehdrste(j,:),'color',ercolors(j,:),...
%         'linewidth',2);
    SLerrorbar(d.time,d.ehdr(j,:),'yError',d.ehdrste(j,:),...
        ['Color=[',num2str(ercolors(j,:)),']'],...
        'linewidth',2,...
        'linesmoothing','on')
end
box off
set(gca,'fontsize',14)
title(['Event-relate plot (',...
    num2str(numVoxels),'voxels ',...
    myROIname,...
    ', r2>',num2str(myr2Thresh),...
    ', hdrlen=',num2str(hdrlen),'vols)',])
ylabel('Signal change (%)')
xlabel('Time (seconds)')

%display everything altogether
img=refreshMLRDisplay(viewGet(getMLRView,'viewNum'));


%number of voxels with r2 greater than different r2 thresholds
%-------------------------------------------------------------
%plot the number of voxels with r2 greater than different thresholds
%change the threshold
myr2Thresholds=0:0.05:0.8;

%count the number of voxels
numVoxelsPerr2Thresh=nan(length(myr2Thresholds),1);
for j=1:length(myr2Thresholds)
    numVoxelsPerr2Thresh(j)=sum(myROI.r2>myr2Thresholds(j));
end

%draw the number of voxels with r2 greater than different thresholds
figure
SLdrawBar(numVoxelsPerr2Thresh,myr2Thresholds,...
    1:length(myr2Thresholds))
xlabel('r2 thresholds')
ylabel('Number of voxels with r2 > threshold')
title([num2str(numVoxelsTotal),'voxels, ',myROIname])


%%draw event-related plot for each voxel
%--------------------------------------

%sort voxels from greatest r2 to lowest
indx=1:numel(myROI.r2);
[~,indx]=sortrows(myROI.r2');
indx=indx(end:-1:1);

%sort time Series
myROI.tSeriesr2Sorted=myROI.tSeries(indx,:);

%figure
myscreensz=get(0,'screensize');
figure('color','w','position',[0.25 0.25 0.25 1].*myscreensz([3 4 3 4]));

for thisVox=1:numVoxelsTotal
    
    %plot this voxel
    subplot( ceil(sqrt(numVoxelsTotal)), ceil(sqrt(numVoxelsTotal)),thisVox)
    
    %deconvolve the hdrs from this voxel's time series
    tSeriesThisVox=myROI.tSeriesr2Sorted(thisVox,:);
    dThisVox=getr2timecourse(tSeriesThisVox,nhdr,hdrlen,scm,viewGet(getMLRView,'framePeriod'));

    %for each stimulus condition (e.g., coherence)
    for j=1:nhdr
        SLerrorbar(dThisVox.time,dThisVox.ehdr(j,:),'yError',dThisVox.ehdrste(j,:),...
            ['Color=[',num2str(ercolors(j,:)),']'],...
            'linewidth',2,...
            'linesmoothing','on',...
            'MarkerSize',5)
        xlim([0 max(dThisVox.time)])
        ylim([min(dThisVox.ehdr(:)) max(dThisVox.ehdr(:))])
    end 
    %axis square
    title(['r2:',num2str(myROI.r2(indx(thisVox)))])
end
















