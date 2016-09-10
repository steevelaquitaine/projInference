%SLfMRIclassificationAnalysis.m
%
%    Author: steeve laquitaine
%      date: 140429 - last modification 150814
%   purpose: fMRI decoding with Fisher linear discriminant
%    
%Description for the GUI:
%
%  >>  mrLoadRet
%
% - (mrLoadRet) File - Analysis - load: /Concatenated/event-related analysis on 
%   concatenated scans 
%
% - (mrLoadRet) set Group: Concatenation
%
% - (mrLoadRet) File - Analysis - Load Concatenated/erAnal/erAnalDir.mat
%   r2 maps appears, now choose variable e.g., "myRandomDir"
%   see http://gru.brain.riken.jp/doku.php/mgl/mglsinglepage?s[]=getstimvol#getstimvol
%
% - load MT's time series ("loadROITSeries.m")
%   tSeries show e.g., 2304 cols (8 scans x (38-2) junk trials x 8 vol.) 
%
% - Get all the voxels (we could also choose the 100 most activated)
%
% - snip up after stimulus and average over stimulus duration (over vols 
%   2:9, i.e., 3 - 13.5s)
%     
% - make one matrix array per ROI (cell array if many ROIs)
%
% - get instances from the ROI (e.g., 144 trials of responses for each category 
%   i{1} or {2} and each trial is a 22 dimensional vector: 22 voxels). i
%   contains the responses.
%
% - leave one out cross-validation for Fisher linear discriminant (default)
%   terminal outputs: "Fisher Classifier produced 78.82% correct and took ....
%   to understand outputs, see: 
%   http://gru.brain.riken.jp/doku.php/mrtools/classificationtools?s[]=leaveoneout


%---------------------
%Initialize parameters
%---------------------
fprintf('%s \n','(SLfMRIclassificationAnalysis) Please set dataPath ...')
dataPath = uigetdir(cd,'Pickup your project e.g., ~/steeve/fMRI_data/s0252015.....');
cd(dataPath)
mrLoadRet([])

%choose "class" (what we decode)
          o.decodedClass = 'myRandomDir';

%choose ROI (where we decode)
             o.myROIname = {'rMT','lMT'}; %or 'V1toMT';
           o.ROIFilePath = '/Users/steeve/data/mlrAnatDB/s0025/mlrROIs';

%choose target task scans
               o.taskNum = 2;
              o.phaseNum = 1;
            o.segmentNum = 2;
o.scanNumFromDescription = 'main task';

%choose flat maps to display on
o.anatFileName = 's0025_right_GM';
o.anatFilePath = '/Users/steeve/data/mlrAnatDB/s0025/surfaces';

%Set lag in vols and # of vols to average together to 
%create BOLD instances. Instances are created by averaging 
%the "blockLen" successive vols starting from "startlag" after
%motion onset
o.startLag = 'startLag=6'; %include 2nd vols and blockLen-1 next vols   
o.blockLen = 'blockLen=24'; %# of averaged vols


%-----------------
%concatenate scans
%-----------------
%check if concatenation group exist
o.groupNames = viewGet(getMLRView,'groupNames');

%if not, concatenate
if sum(strcmp(groupNames,'Concatenation'))~=1;

    %set motion compensation group for concatenation
    curGroup = 'MotionComp';
    viewSet(getMLRView,'curGroup',curGroup); %set motion comp group
    scanNum = viewGet(getMLRView,'scanNumFromDescription',scanNumFromDescription,curGroup);
    
    %detrend and high-pass filter the scans you select, get rid of any
    %junk frames, do a de-trending step using a hi-pass filter if you
    %specify that, and then concatenate the resulting time series into a
    %single time series that will be saved into a new group called
    %"Concatenation". The .mat files are saved in folder "Concatenation".
    [~,params] = concatTSeries(getMLRView,[],'justGetParams=1',... %set concatenation parameters
        'defaultParams=1','scanList',scanNum);
    
    %concatenate (stack preprocessed scans into one long scan)
    concatTSeries(getMLRView,params);
    fprintf('%s \n','(SLfMRIclassificationAnalysis) Motion compensation scans have been concatenated.')    
else
    fprintf('%s \n','(SLfMRIclassificationAnalysis) Motion compensation scans were already concatenated.')
end


%choose concatenated group and scan for classification analysis
viewSet(getMLRView,'curGroup','Concatenation');
nScans = viewGet(getMLRView,'nScans');
viewSet(getMLRView,'curScan',nScans);

%---------------------
%Task conditions found
%---------------------

%HERE !!!!!!!!!!!!
%Select vols for selected task conditions (we can input the task conditions wanted)


%---------
% CLASSIFY
%---------

%load flat map anatomy (requires 2 files e.g., 's025_leftFlat.hdr/img/mat')
%some flat maps have been created for retinotopy. We can just use those
%those too.
loadAnat(getMLRView,anatFileName,anatFilePath);

%get acquired volumes for selected class at motion onset (segment 2)
o.stimvol = getStimvol(getMLRView,decodedClass,'taskNum',taskNum,'phaseNum',phaseNum,'segmentNum',segmentNum);

%load ROI
loadROI(getMLRView,myROIname,0,ROIFilePath)

%display in mrLoadRet GUI 
img = refreshMLRDisplay(viewGet(getMLRView,'viewNum'));

%load ROI time series
o.myROI = loadROITSeries(getMLRView,myROIname);

%---------------------------------------------
%Accuracy for full period after stimulus onset
%---------------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset.
o.myROI = slGetInstances(getMLRView,o.myROI,o.stimvol,'n=inf',o.startLag,blockLen);
o.myROI = o.myROI{1};
i = myROI.classify.instances;
data.fullPeriod = leaveOneOut(i);
%weightMap=getClassifierWeightMap(getMLRView,myROI,stimvol);
% 
% %remove current analysis and look at the classifier's weightMap
% % - view menu - remove Analysis
% mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
%     'colormapType','normal','cmap',splitcmap,...
%     'range=[-9 9]','clip=[-9 9]');

%CLASSIFY z-scored data
[izscored,pSettings] = preprocessInstances(i,'zscore=1');
data.Zscored.fullPeriod = leaveOneOut(izscored);

%--------------------------------------------------------------
%Visualize accuracy for different periods after stimulus onset
%--------------------------------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset. Set lag in volumes and length of block data  in vol. over which to 
%average data into instances.
startLags = {'startLag=1','startLag=2','startLag=3','startLag=4','startLag=1',...
    'startLag=2','startLag=3','startLag=4','startLag=1'};
blockLens = {'blockLen=2','blockLen=2','blockLen=2','blockLen=2','blockLen=6',...
    'blockLen=6','blockLen=6','blockLen=6','blockLen=10'};
tic
for j = 1 : numel(startLags)
    
    %set instance averaging period
    startLag = startLags{j};
    blockLen = blockLens{j};
    
    %get instances
    myROI = slGetInstances(getMLRView,myROI,stimvol,'n=inf',startLag,blockLen);
    myROI = myROI{1};
    i = myROI.classify.instances;
    
    %classify
    data.perperiods{j} = leaveOneOut(i);
    %weightMap=getClassifierWeightMap(getMLRView,myROI,stimvol);
    %
    % %remove current analysis and look at the classifier's weightMap
    % % - view menu - remove Analysis
    % mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
    %     'colormapType','normal','cmap',splitcmap,...
    %     'range=[-9 9]','clip=[-9 9]');
    
    %CLASSIFY z-scored data
    [izscored,pSettings] = preprocessInstances(i,'zscore=1');
    data.Zscored.perperiods{j} = leaveOneOut(izscored);
    
    %clear up stuffs
    myROI.classify = [];
    
    %print
    fprintf('%12s /n',['(SLfMRIclassificationAnalysis) ',num2str(round(j/numel(startLags)*100)),'% complete'])
end
toc


%------
%backup
%------
mfilename = SLgetActivemFile;
myMfile = SLbackupMfileInWSpace(mfilename);
o.mfilename  = mfilename;
o.myMfile = myMfile;

%backup m.file with date and time
filename = matlab.desktop.editor.getActiveFilename;
SLBackup(filename)

