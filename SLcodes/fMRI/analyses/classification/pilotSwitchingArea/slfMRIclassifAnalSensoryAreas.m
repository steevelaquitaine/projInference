
%slfMRIclassifAnalSensoryAreas.m
%
%    Author: steeve laquitaine
%      date: 150922
%   purpose: decoding motion direction or coherence from ROI voxels Bold with Fisher linear discriminant 
%    
%     usage: 
%                o.myPath = '~/dataold/datafMRI/sltaskdotdirfmri05/s02520150814'
%               o.myGroup = 'Concatenation';
%        o.myGroupScanNum = 1;
%              o.scanList = 1:18;                %for concatenation if needed
%               o.myClass = 'myRandomDir';              %only 6% coh
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
%
%note:
%
%   Duration : about 36 min for 18 x 614 vols.

function o = slfMRIclassifAnalSensoryAreas(o)

tic 

%launch mrLoadRet
cd(o.myPath)
v = mrLoadRet;

%Sum info
fprintf('----------------------Data info ----------------------------------- \n')
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) ROI: ',o.myROIname)
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task: ',o.taskNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task segment: ',o.segmentNum)
fprintf('%s %s %s %s %s %s \n \n','(slfMRIclassifAnalSensoryAreas) Vols: ',o.startLag,' to ',o.startLag,'+',o.blockLen)

fprintf('----------------------Classification info------------------------- \n')
fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Fisher linear discriminant')
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Class:',o.myClass)
fprintf('------------------------------------------------------------------ \n')

%choose concatenated group and scan
viewSet(getMLRView,'curGroup',o.myGroup);

%-----------------
%CLASSIFY from ROI
%-----------------
%Vols used to decode class at event (Nclasses cells of 1 by Ntrialrep vol #)
[o.stimvol o.taskCond] = getStimvol(v,o.myClass,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
fprintf('%s \n',' (slfMRIclassifAnalSensoryAreas) Task conditions found are: ')
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

keyboard

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
fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Now decoding at contiguous time steps....')

%-----------------------------------------
%Plot accuracy by period after event onset
%-----------------------------------------
%Create instances by averaging volumes in a defined period after stimulus 
%onset. Set lag in volumes and length of block data  in vol. over which to 
%average data into instances.
%cover a up to 20 sec after event onset with time steps of 2 vols. 
o.tr = 0.5;
for i = 1 : round(20/o.tr)
    startLags{i} = ['startLag=' num2str(i)];
end
blockLens   = 'blockLen=2';
myROI       = o.myROI;
stimvol     = o.stimvol;
lengthsteps = numel(startLags);

parfor j = 1 : numel(startLags)
    
    %classify raw instances
    %----------------------
    fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Decoding Raw Bold....')
    myROIj{j}     = getInstances(v,myROI,stimvol,'n=inf',startLags{j},blockLens);
    ij{j}         = myROIj{j}{1}.classify.instances;
    perperiods{j} = leaveOneOut(ij{j});
    
    %plot MAP of voxels weight
    %o.weightMap = getClassifierWeightMap(v,o.myROI,stimvol);
    
    %remove current analysis and look at the classifier's weightMap
    % % - view menu - remove Analysis
    %mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
    %     'colormapType','normal','cmap',splitcmap,...
    %     'range=[-9 9]','clip=[-9 9]');
    
    %classify z-scored instances
    %---------------------------
    fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Decoding z-scored Bold....')
    [izscoredj{j},pSettingsj{j}] = preprocessInstances(ij{j},'zscore=1');
    ZscoredPerperiods{j}         = leaveOneOut(izscoredj{j});
    
    %print
    fprintf('%12s %i %s \n','(slfMRIclassifAnalSensoryAreas) ',round(j/lengthsteps *100),'% complete')
    myROIj = [];
end

%output
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

%plot accuracy from raw instances
SLerrorbar([1:lengthsteps]*o.tr, o.classifRes.raw.perperiods.Accuracy, 'yError',o.classifRes.raw.perperiods.AccuracySTE,...
['Color=[',num2str([.8 .8 .8]),']'],'MarkerSize=7','linewidth',2,'linesmoothing','on');

%plot accuracy from zscored instances
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

fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Accuracy Raw data     : ',o.classifRes.raw.fullPeriod.correct)
fprintf('%s %i \n \n','(slfMRIclassifAnalSensoryAreas) Accuracy zscored data : ',o.classifRes.Zscored.fullPeriod.correct)



fprintf('----------------------Results classif (per period) ------------------------- \n')

fprintf('%s %.2f \n','(slfMRIclassifAnalSensoryAreas) Accuracy Raw data     : ',o.classifRes.raw.perperiods.Accuracy)
fprintf('%s %.2f \n \n','(slfMRIclassifAnalSensoryAreas) Accuracy zscored data : ',o.classifRes.Zscored.perperiods.Accuracy)

%backup
mkdir(o.ROIname)
cd(o.ROIname)
save(['Classif' o.myClass o.ROIname date],'o')
saveas(gcf, ['Classif' o.myClass o.ROIname], 'fig')
o.duration = toc;
