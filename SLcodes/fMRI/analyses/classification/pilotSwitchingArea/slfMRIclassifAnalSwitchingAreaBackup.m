
%slfMRIclassifAnalSwitchingArea.m
%
%    Author: steeve laquitaine
%      date: 150922
%   purpose: decoding switching (1 or 0 | near prior or like) from ROI voxels Bold with Fisher linear discriminant
%    
%     usage: 
%                o.myPath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814'
%               o.myGroup = 'MotionComp';%'Concatenation'
%        o.myGroupScanNum = 1;
%              o.scanList = 1:18;                                %for concatenation if needed
%               o.myClass = 'mySwitch';                          %only 6% coh
%             o.ClassType = 'svm';  
%             o.myROIname = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrROIs/V1.mat';%Currently selected ROI (GUI)
%               o.taskNum = 2;%choose target task scans
%              o.phaseNum = 1;
%            o.segmentNum = 2;
%o.scanNumFromDescription = 'main task';
%          o.anatFileName = 's0025_flatL_WM_occipital_Rad90';%choose flat maps to display on
%          o.anatFilePath = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
%              o.startLag = 'startLag=1'; %include 2nd vols and blockLen-1 next vols   
%              o.blockLen = 'blockLen=10'; %# of averaged vols

function o = slfMRIclassifAnalSwitchingArea(o)

tic

%launch mrLoadRet
cd(o.myPath)
v = mrLoadRet

%Sum info
fprintf('----------------------Data info ----------------------------------- \n')
fprintf('%s %s \n','(SLfMRIclassifAnalSwitchingArea) ROI: ',o.myROIname)
fprintf('%s %s \n','(SLfMRIclassifAnalSwitchingArea) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(SLfMRIclassifAnalSwitchingArea) Task: ',o.taskNum)
fprintf('%s %i \n','(SLfMRIclassifAnalSwitchingArea) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(SLfMRIclassifAnalSwitchingArea) Task segment: ',o.segmentNum)
fprintf('%s %s %s %s %s %s \n \n','(SLfMRIclassifAnalSwitchingArea) Vols: ',o.startLag,' to ',o.startLag,'+',o.blockLen)

fprintf('----------------------Classification info------------------------- \n')
fprintf('%s \n','(SLfMRIclassifAnalSwitchingArea) Fisher linear discriminant')
fprintf('%s %s \n','(SLfMRIclassifAnalSwitchingArea) Class:',o.myClass)
fprintf('------------------------------------------------------------------ \n')

%choose concatenated group and scan
viewSet(getMLRView,'curGroup',o.myGroup);

keyboard


%%%--------------------------------------
%Create "Switching" classification variable
%%%--------------------------------------
%load beh. data
data = loadScan(v);
e    = getTaskParameters(data.stimfile{1}.myscreen,data.stimfile{1}.task);
e    = e{o.taskNum};  
[~,data.estimatesDeg] = SLcart2polar(cell2mat(e.randVars.prodcoor')); %estimates

%create variable indicating when subject responded
missing         = find(isnan(data.estimatesDeg)); 
isresponse      = ones(length(data.estimatesDeg),1);
isresponse(missing) = 0;

%label trials as "switched-to-prior" (1) and "switched-to-motion dir" (2)
data.motiondir  = e.randVars.myRandomDir;  					          %direction
distPandL = abs([SLvectors2signedAngle(data.estimatesDeg,225,'polar') SLvectors2signedAngle(data.estimatesDeg,data.motiondir,'polar')])
mySwitch = nan(length(data.estimatesDeg),1)
for i = 1 : length(data.estimatesDeg)
	[~,mySwitch(i,1)] = min(distPandL(i,:)); %label switch 
end
mySwitch(missing) = nan; 		%label missing
nMiss = sum(isnan(mySwitch));
nTrials  = length(mySwitch);
fprintf('%s %i %s %i %s %i %s \n','(SLfMRIclassifAnalSwitchingArea) ',nMiss/nTrials*100,'% (' ,nMiss,'/',nTrials,') of the responses are missing.')

%if raw Stim has not been backed up yet create a backup
cd(o.myPath)
cd Etc
[pthi stimFilenm] = fileparts(data.stimfile{1}.filename);
if ~exist([pthi '/rawStimFiles/' stimFilenm '.mat'],'file')
	mkdir rawStimFiles
	copyfile(data.stimfile{1}.filename,'rawStimFiles/') %backup raw file
    cd ..
    fprintf('%s %s %s \n','(SLfMRIclassifAnalSwitchingArea)',stimFilenm, 'has been backed up in rawStimFiles.')    
end

%------------------------------
%ADD NEW VARIABLES TO STIM FILE
%------------------------------
%add variables to new stimfile (for fMRI analyses) ; create folder and backup raw stimfiles !!!
%add "isresponse" and "switch" variables
load(stimFilenm)
task{o.taskNum}{1}.randVars.isresponse = isresponse;                  
task{o.taskNum}{1}.randVars.mySwitch   = mySwitch;                   
nAddedVar   = 2;
addedVarnm  = {'isresponse','mySwitch'};
addedVarLen = [nTrials nTrials];
for i = 1 : nAddedVar
    task{o.taskNum}{1}.randVars.names_(task{o.taskNum}{1}.randVars.n_ + i) = addedVarnm(i);
    task{o.taskNum}{1}.randVars.varlen_(task{o.taskNum}{1}.randVars.n_ + i) = addedVarLen(i);
end
task{o.taskNum}{1}.randVars.n_   = task{o.taskNum}{1}.randVars.n_ + nAddedVar;
save(stimFilenm,'fixStimulus','myscreen','stimulus','task')                   
load(stimFilenm)

%-----------------
%CLASSIFY from ROI
%-----------------
%Vols used to decode class at event (Nclasses cells of 1 by Ntrialrep vol #)
[o.stimvol o.taskCond] = getStimvol(v,{{'mySwitch=1'},{'mySwitch=2'}},'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
fprintf('%s \n',' (SLfMRIclassifAnalSwitchingArea) Task conditions found are: ')
for i = 1 : length(o.taskCond)
    disp(o.taskCond{i})
end
o.myROI = loadROITSeries(v,o.myROIname);

%Accuracy for full stim segment
%------------------------------
%Create instances by averaging vols in a period after stim. onset
o.myROI = getInstances(v,o.myROI,o.stimvol,'n=inf',o.startLag,o.blockLen);
o.myROI = o.myROI{1};
o.i     = o.myROI.classify.instances; %Bold average over vols


%Fisher linear discriminant
%--------------------------

%RAW data
%--------
o.classifRes.raw.fullPeriod = leaveOneOut(o.i);

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
o.classifRes.Zscored.fullPeriod = leaveOneOut(o.izscored);

%--------------------------------------------------------------
%Visualize accuracy for different periods after stimulus onset
%--------------------------------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset. Set lag in volumes and length of block data  in vol. over which to 
%average data into instances.
fprintf('%s \n','(SLfMRIclassifAnalSwitchingArea) Now decoding at contiguous time steps....')

%--------------------------------------------------------------
%Visualize accuracy for different periods after stimulus onset
%--------------------------------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset. Set lag in volumes and length of block data  in vol. over which to 
%average data into instances.
%o.startLags = {'startLag=1','startLag=2','startLag=3','startLag=4','startLag=1',...
%    'startLag=2','startLag=3','startLag=4'};

%o.blockLens = {'blockLen=2','blockLen=2','blockLen=2','blockLen=2','blockLen=6',...
%    'blockLen=6','blockLen=6','blockLen=6'};

%cover a up to 20 sec after event onset with time steps of 2 vols 
o.tr = 0.5;
for i = 1 : round(20/o.tr)
    startLags{i} = ['startLag=' num2str(i)];
end
blockLens = {'blockLen=2'};

myROI   = o.myROI;
stimvol = o.stimvol;
lengthsteps = numel(startLags);

parfor j = 1 : numel(startLags)
    
    %Instances
    %---------
    startLag = startLags{j}; % start vols
    blockLen = blockLens; % # of averaged vols
    myROIj   = getInstances(v,myROI,stimvol,'n=inf',startLag,blockLen);
    ij{j}    = myROIj{1}.classify.instances;
    
    %Classify
    %--------
    fprintf('%s \n','(SLfMRIclassifAnalSwitchingArea) Decoding Raw Bold....')
    perperiods{j} = leaveOneOut(ij{j});
    
    %plot MAP of voxels weight
    %o.weightMap = getClassifierWeightMap(v,o.myROI,stimvol);
    
    %remove current analysis and look at the classifier's weightMap
    % % - view menu - remove Analysis
    %mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
    %     'colormapType','normal','cmap',splitcmap,...
    %     'range=[-9 9]','clip=[-9 9]');
    
    %CLASSIFY z-scored data
    fprintf('%s \n','(SLfMRIclassifAnalSwitchingArea) Decoding z-scored Bold....')
    [izscoredj,pSettingsj] = preprocessInstances(ij{j},'zscore=1');
    ZscoredPerperiods{j} = leaveOneOut(izscoredj);
        
    myROIj = [];

    %print
    fprintf('%12s %i %s \n','(SLfMRIclassifAnalSwitchingArea) ',round(j/lengthsteps *100),'% complete')
end

o.startLags = startLags;
o.blockLens = blockLens;
o.ij        = ij;

%---------------------------
%Plot accuracy by time steps
%---------------------------
figure('color','w')
%get accuracies
for i = 1:lengthsteps
    o.classifRes.raw.perperiods.Accuracy(i)     = perperiods{i}.correct;
    o.classifRes.raw.perperiods.AccuracySTE(i)  = perperiods{i}.correctSTE;
    o.classifRes.Zscored.perperiods.Accuracy(i) = ZscoredPerperiods{i}.correct;
    o.classifRes.Zscored.perperiods.AccuracySTE(i) = ZscoredPerperiods{i}.correctSTE;
end

SLerrorbar([1:lengthsteps]*o.tr, o.classifRes.raw.perperiods.Accuracy, 'yError',o.classifRes.raw.perperiods.AccuracySTE,...
['Color=[',num2str([.8 .8 .8]),']'],'MarkerSize=7','linewidth',2,'linesmoothing','on');

hold on; SLerrorbar([1:lengthsteps]*o.tr, o.classifRes.Zscored.perperiods.Accuracy, 'yError',o.classifRes.Zscored.perperiods.AccuracySTE,...
['Color=[',num2str([0 0 0]),']'],'MarkerSize=7','linewidth',2,'linesmoothing','on');

legend('raw Bold (grey)','z-scored Bold (black)')

%full backup
o.perperiods{i}.correct = perperiods{i}.correct;
o.perperiods{i}.correctSTE = perperiods{i}.correctSTE;
o.ZscoredPerperiods{i}.correct = ZscoredPerperiods{i}.correct;
o.ZscoredPerperiods{i}.correctSTE = ZscoredPerperiods{i}.correctSTE;

%get task segments and add them to plot
a = viewGet(v,'stimfile')
for i = 1: length(a)
    seglens(i,:) = mean([a{i}.task{o.taskNum}{1}.segmin; a{i}.task{o.taskNum}{1}.segmax]); %average segment lengths
end
seglens = unique(seglens,'rows');
if size(seglens,1) ~= 1
    fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) The scans trials have different segment lengths. e.g., trial 1 : [1 0.3 5 3 11.7] and trial N : [1 3 10 3 20]. The code currently requires all scans to have same segment lengths across trials. e.g., all trials : [1 0.3 5 3 11.7]')
    keyboard
end
cumsumseglen = [0 cumsum(seglens)];
TrialEventTimes = cumsumseglen - cumsumseglen(o.segmentNum);
xlim([min(TrialEventTimes)-1 max([TrialEventTimes lengthsteps*o.tr])])
eventcolors = linspecer(length(seglens));
for i = 1 : length(seglens)+1
    vline(TrialEventTimes(i))
end

% with info
o.nClasses = length(stimvol);
randAccuracy = 1./o.nClasses;
hline(randAccuracy)
nVoxels = size(o.i{1},2);
for i = 1 : o.nClasses
     o.repeatPerClass(i) = size(stimvol{i},2);
end
[a o.ROIname]=fileparts(o.myROIname);
title([o.ROIname ' - ' o.myClass ' (' num2str(o.nClasses) ' classes) - ' num2str(nVoxels) ' voxels - ' num2str(min(o.repeatPerClass)) ' reps min per class'])
ylim([randAccuracy+[ -0.2 .2]])
xlabel('Time after event (sec)')
ylabel('Accuracy (% correct)')
box off

%-------
%Summary
%-------
fprintf('----------------------Results classif. (full period) ----------------------- \n')

fprintf('%s %i \n','(SLfMRIclassifAnalSwitchingArea) Accuracy Raw data     : ',o.classifRes.raw.fullPeriod.correct)
fprintf('%s %i \n \n','(SLfMRIclassifAnalSwitchingArea) Accuracy zscored data : ',o.classifRes.Zscored.fullPeriod.correct)



fprintf('----------------------Results classif (per period) ------------------------- \n')

fprintf('%s %.2f \n','(SLfMRIclassifAnalSwitchingArea) Accuracy Raw data     : ',o.classifRes.raw.perperiods.Accuracy)
fprintf('%s %.2f \n \n','(SLfMRIclassifAnalSwitchingArea) Accuracy zscored data : ',o.classifRes.Zscored.perperiods.Accuracy)

o.duration = toc;

%backup
mkdir(o.ROIname)
cd(o.ROIname)
save(['ClassifSwitch' o.myClass o.ROIname date],'o')
saveas(gcf, ['ClassifSwitch' o.myClass o.ROIname], 'fig')
