%SLfMRIclassifAnalSensAreaGui.m
%
%    Author: steeve laquitaine
%      date: 140429 - last modification 150814
%   purpose: fMRI decoding of physical motion stimulus direction with Fisher linear discriminant
%    
%     usage: 
%               1) mrLoadRet
%               2) choose your stuff (group, ROI, Scan, Base)
%               3) Plots - Interrogator
%               4) Write down : SLfMRIclassifAnalSensAreaGui in Interrogator box
%               5) Click on you ROI
%
%            
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

%o = SLfMRIclassifAnalSensAreaGui(view,overlayNum,scan,x,y,s,roi,varargin)

function o = SLfMRIclassifAnalSensAreaGui(view,overlayNum,scan,x,y,s,roi)

% check arguments
%if ~any(nargin == [0])
%  help SLfMRIclassifAnalSensAreaGui
%  return
%end

v = getMLRView;

%---------------------
%Initialize parameters
%---------------------
%fprintf('%s \n','(SLfMRIclassificationAnalysis) Please set dataPath ...')
%dataPath = uigetdir(cd,'Pickup your project e.g., ~/steeve/fMRI_data/s0252015.....');
%cd(dataPath)
%mrLoadRet([])

%choose "class" (what we decode)
          %o.myClass = 'myRandomDir';
          %o.myClass = 'myRandomDir_x_myRandomCoh';
             o.myClass = 'myRandomDir_x_myRandomCoh=0.12';       %only 6% coh
          %o.myClass = {{'myRandomDir_x_myRandomCoh=0.06'},{'myRandomDir_x_myRandomCoh=0.12'}}

%choose classifier
            o.ClassType = 'svm';  

%Currently selected ROI (GUI)
             o.myROIname = v.ROIs(v.curROI).name;

           %o.ROIFilePath = '/Users/steeve/data/mlrAnatDB/s0025/mlrROIs';

%choose target task scans
               o.taskNum = 2;
              o.phaseNum = 1;
            o.segmentNum = 2;
o.scanNumFromDescription = 'main task';

%choose flat maps to display on
%o.anatFileName = 's0025_right_GM';
o.anatFileName = v.baseVolumes.name;
%o.anatFilePath = '/Users/steeve/data/mlrAnatDB/s0025/surfaces';

%Set lag in vols and # of vols to average together to 
%create BOLD instances. Instances are created by averaging 
%the "blockLen" successive vols starting from "startlag" after
%motion onset
o.startLag = 'startLag=1'; %include 2nd vols and blockLen-1 next vols   
o.blockLen = 'blockLen=10'; %# of averaged vols

%Sum info
fprintf('----------------------Data info ----------------------------------- \n')
fprintf('%s %s \n','(SLfMRIclassifAnalSensAreaGui) ROI: ',o.myROIname)
fprintf('%s %s \n','(SLfMRIclassifAnalSensAreaGui) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(SLfMRIclassifAnalSensAreaGui) Task: ',o.taskNum)
fprintf('%s %i \n','(SLfMRIclassifAnalSensAreaGui) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(SLfMRIclassifAnalSensAreaGui) Task segment: ',o.segmentNum)
fprintf('%s %s %s %s %s %s \n \n','(SLfMRIclassifAnalSensAreaGui) Vols: ',o.startLag,' to ',o.startLag,'+',o.blockLen)

fprintf('----------------------Classification info------------------------- \n')
fprintf('%s \n','(SLfMRIclassifAnalSensAreaGui) Fisher linear discriminant')
fprintf('%s %s \n','(SLfMRIclassifAnalSensAreaGui) Class:',o.myClass)
fprintf('------------------------------------------------------------------ \n')


%-----------------
%concatenate scans
%-----------------
%check if concatenation group exist
%o.groupNames = viewGet(getMLRView,'groupNames');

%if not, concatenate
%if sum(strcmp(groupNames,'Concatenation'))~=1;

    %set motion compensation group for concatenation
%    curGroup = 'MotionComp';
%    viewSet(getMLRView,'curGroup',curGroup); %set motion comp group
%    scanNum = viewGet(getMLRView,'scanNumFromDescription',scanNumFromDescription,curGroup);
    
    %detrend and high-pass filter the scans you select, get rid of any
    %junk frames, do a de-trending step using a hi-pass filter if you
    %specify that, and then concatenate the resulting time series into a
    %single time series that will be saved into a new group called
    %"Concatenation". The .mat files are saved in folder "Concatenation".
%    [~,params] = concatTSeries(getMLRView,[],'justGetParams=1',... %set concatenation parameters
%        'defaultParams=1','scanList',scanNum);
    
    %concatenate (stack preprocessed scans into one long scan)
%    concatTSeries(getMLRView,params);
%    fprintf('%s \n','(SLfMRIclassificationAnalysis) Motion compensation scans have been concatenated.')    
%else
%    fprintf('%s \n','(SLfMRIclassificationAnalysis) Motion compensation scans were already concatenated.')
%end


%choose concatenated group and scan for classification analysis
%viewSet(getMLRView,'curGroup','Concatenation');
o.nScans = viewGet(v,'nScans');

%---------
% CLASSIFY
%---------
%Vols used to decode class at motion onset (Nclasses cells of 1 by Ntrialrep vol #)
[o.stimvol o.taskCond] = getStimvol(v,o.myClass,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);

fprintf('%s \n',' (SLfMRIclassifAnalSensAreaGui) Task conditions found are: ')
for i = 1 : length(o.taskCond)
    disp(o.taskCond{i})
end

%ROI time series
o.myROI = loadROITSeries(v,o.myROIname);

%Accuracy for full stim segment
%------------------------------
%Create instances by averaging vols in a period after stim. onset
o.myROI = getInstances(v,o.myROI,o.stimvol,'n=inf',o.startLag,o.blockLen);
o.myROI = o.myROI{1};
o.i = o.myROI.classify.instances; %Bold average over vols


%Fisher linear discriminant
%--------------------------

%RAW data
%--------
o.data.fullPeriod = leaveOneOut(o.i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO DO
% - Test other classifiers HERE
%       - SVM
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Map of voxels weight in classification
%--------------------------------------
%o.weightMap = getClassifierWeightMap(v,o.myROI,o.stimvol);
%mrDispOverlay(o.weightMap,o.myROI.scanNum,o.myROI.o.groupNum,v,...
%     'colormapType','normal','cmap',splitcmap,...
%     'range=[-9 9]','clip=[-9 9]');

%z-scored data
%-------------
[o.izscored,o.pSettings] = preprocessInstances(o.i,'zscore=1');
o.data.Zscored.fullPeriod = leaveOneOut(o.izscored);

%--------------------------------------------------------------
%Visualize accuracy for different periods after stimulus onset
%--------------------------------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset. Set lag in volumes and length of block data  in vol. over which to 
%average data into instances.

o.startLags = {'startLag=1','startLag=2','startLag=3','startLag=4','startLag=1',...
    'startLag=2','startLag=3','startLag=4'};

o.blockLens = {'blockLen=2','blockLen=2','blockLen=2','blockLen=2','blockLen=6',...
    'blockLen=6','blockLen=6','blockLen=6'};

for j = 1 : numel(o.startLags)
    
    %Instances
    %---------
    startLag = o.startLags{j}; % start vols
    blockLen = o.blockLens{j}; % # of averaged vols
    o.myROIj = slGetInstances(v,o.myROI,o.stimvol,'n=inf',startLag,blockLen);
    o.ij     = o.myROIj{1}.classify.instances;
    
    %Classify
    %--------
    o.data.perperiods{j} = leaveOneOut(o.ij);
    
    %plot MAP of voxels weight
    %o.weightMap = getClassifierWeightMap(v,o.myROI,stimvol);
    
    %remove current analysis and look at the classifier's weightMap
    % % - view menu - remove Analysis
    %mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
    %     'colormapType','normal','cmap',splitcmap,...
    %     'range=[-9 9]','clip=[-9 9]');
    
    %CLASSIFY z-scored data
    [o.izscoredj,o.pSettingsj] = preprocessInstances(o.ij,'zscore=1');
    data.Zscored.perperiods{j} = leaveOneOut(o.izscoredj);
        
    o.myROIj = [];

    %print
    fprintf('%12s %i %s \n','(SLfMRIclassifAnalSensAreaGui) ',round(j/numel(o.startLags)*100),'% complete')
end

%Sum Results
fprintf('----------------------Results classif. (full period) ----------------------- \n')
fprintf('%s %i \n','(SLfMRIclassifAnalSensAreaGui) Accuracy Raw data     : ',o.data.fullPeriod.correct)
fprintf('%s %i \n \n','(SLfMRIclassifAnalSensAreaGui) Accuracy zscored data : ',o.data.Zscored.fullPeriod.correct)

fprintf('----------------------Results classif (per period) ------------------------- \n')
fprintf('%s \n','(SLfMRIclassifAnalSensAreaGui) Accuracy')
fprintf('------------------------------------------------------------------------------ \n')

save('ClassifResults','o')
