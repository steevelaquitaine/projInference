
%slfMRIclassifAnalSwitchingArea2.m
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

function [o,c,v] = slfMRIclassifAnalSwitchingArea2(o,varargin)

%pick time        
if sum(strcmp(varargin,'accuracyByTime'))
    
    %accuracy time course
    o.myAnalysis = 'accByTime';
        
elseif sum(strcmp(varargin,'accuracyAtTime'))
    
    %accuracy at specific time
    o.myAnalysis = 'accAtTime';
    pos = find(strcmp(varargin,'accuracyAtTime')) + 1;
    %[start end] vols to average for classification
    o.volsToClass = varargin{pos};        
    
end

%launch clean mrLoadRet 
tic
cd(o.myPath)
v = mrLoadRet([]);

%--------
%Sum info
%--------
fprintf('----------------------Data info ----------------------------------- \n')
nRois = length(o.myROIname);
for i = 1 : nRois
    fprintf('%s %s \n','(slfMRIclassifAnalSwitchingArea2) ROI: ',o.myROIname{i})
end
fprintf('%s %s \n','(slfMRIclassifAnalSwitchingArea2) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(slfMRIclassifAnalSwitchingArea2) Task: ',o.taskNum)
fprintf('%s %i \n','(slfMRIclassifAnalSwitchingArea2) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(slfMRIclassifAnalSwitchingArea2) Task segment: ',o.segmentNum)

fprintf('----------------------Classification info------------------------- \n')
fprintf('%s \n','(slfMRIclassifAnalSwitchingArea2) Fisher linear discriminant')
fprintf('%s %s \n','(slfMRIclassifAnalSwitchingArea2) Class:',o.myClass)
fprintf('------------------------------------------------------------------ \n')

%group & scan
v = viewSet(getMLRView,'curGroup',o.myGroup);
v = viewSet(getMLRView,'curScan',o.myScan);
refreshMLRDisplay(v);                     %update display

%---------------------------------------------------------------
%loop over stimfiles and get data and task conditions.
%Update them with new variable (switching and missing responses)
%---------------------------------------------------------------
myStimfiles = viewGet(v,'stimfile');
nStimf = length(myStimfiles);
fprintf('%s %i %s \n','(slfMRIclassifAnalSwitchingArea2)',nStimf, 'stim files have been found.')    
fprintf('%s \n','%--------------------------------------- ----------------------------------------------------')    
fprintf('%s \n','%--------------(slfMRIclassifAnalSwitchingArea2: add new variables to stimfiles)-----------------')    
fprintf('%s \n','%--------------------------------------- ----------------------------------------------------')    
for stimi = 1 :  nStimf

    %%%--------------------------------------
    %Create "Switching" classification variable
    %%%--------------------------------------
    %load beh. data
    %data = loadScan(v); memory load
    e = getTaskParameters(myStimfiles{stimi}.myscreen,myStimfiles{stimi}.task);
    e = e{o.taskNum};  
    [~,estimatesDeg] = SLcart2polar(cell2mat(e.randVars.prodcoor')); %estimates

    %create variable indicating when subject responded
    missing         = find(isnan(estimatesDeg)); 
    isresponse      = ones(length(estimatesDeg),1);
    isresponse(missing) = 0;

    %label trials as "switched-to-prior" (1) and "switched-to-motion dir" (2)
    motiondir  = e.randVars.myRandomDir;  					          %direction
    distPandL = abs([SLvectors2signedAngle(estimatesDeg,225,'polar') SLvectors2signedAngle(estimatesDeg,motiondir,'polar')]);
    mySwitch = nan(length(estimatesDeg),1);
    for i = 1 : length(estimatesDeg)
    	[~,mySwitch(i,1)] = min(distPandL(i,:)); %label switch 
    end
    mySwitch(missing) = nan; 		%label missing
    nMiss = sum(isnan(mySwitch));
    nTrials  = length(mySwitch);
    fprintf('%s %.2f %s %i %s %i %s \n','(slfMRIclassifAnalSwitchingArea2) ',nMiss/nTrials*100,'% (' ,nMiss,'/',nTrials,') of the responses are missing.')

    %if raw Stim has not been backed up yet create a backup
    cd(o.myPath)
    cd Etc
    [pthi stimFilenm] = fileparts(myStimfiles{stimi}.filename);
    if ~exist([pthi '/rawStimFiles/' stimFilenm '.mat'],'file')
    	mkdir rawStimFiles
    	copyfile(myStimfiles{stimi}.filename,'rawStimFiles/') %backup raw file
        %cd ..
        fprintf('%s %s %s \n','(slfMRIclassifAnalSwitchingArea2)',stimFilenm, 'has been backed up in rawStimFiles.')    
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

    fprintf('%s \n','(slfMRIclassifAnalSwitchingArea2) Adding variables "isresponse" and "mySwitch" to the stimfiles')    
    fprintf('%s \n','%--------------------------------------- done ----------------------------------------------------')    
end


%----------
%INFO CHECK
%----------
%Vols to decode class at event (Nclasses cells of 1 by Ntrialrep vol #)
[o.stimvol, o.taskCond] = getStimvol(v,o.myClass,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
fprintf('%s %i \n','---------------- Please check that those info are correct -----------------------------')
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Current group: ', viewGet(v,'curGroup'))
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Current scan: ', viewGet(v,'curScan'))
fprintf('%s %i %s \n','(slfMRIclassifAnalSensoryAreas) Max # of vol found : ',max([o.stimvol{:}]),'vols')
fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Task conditions found are: ')
for i = 1 : length(o.taskCond)
    disp(o.taskCond{i})
end
fprintf('------------------------------------------------------------- \n')

%get rid of NaN trials
for i = 1 : length(o.stimvol)
    if isempty(o.stimvol{i})
        o.stimvol(i) = [];
        o.taskCond(i) = [];
    end
end

%------------------------- load ROITseries -------------------------------
%This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
%GET rid of it asap. Prob is er analysis doesn't load and store the
%entire volume but load slice by slice. Thus cannot be used here to get
%the ROI TSeries. The fastest way so far is to save the ROI TSeries once
%run once and reload them instead of loadROITSeries them.
%e.myROI.scanCoords = getROICoordinates(v,myROIname);
o.YesorNo = input(['!! README !! Loading ROI tseries takes a while for long'...
    ' scans (e.g., concatenation). Do you want to load existing ROItseries ? [y/n]'],'s');

o.nRois = length(o.myROIname);

%-----------------
%CLASSIFY from ROI
%-----------------
%repeat analyses by rois
for ROIi = 1 : nRois
    
    %first cleanup everything
    keep = {'v','o','ROIi','keep','a','i'}; a = whos;
    for i = 1 : length(a)
        if ~sum(strcmp(a(i).name,keep))>0
            clear(a(i).name)
        end
    end
    clear keep
    clear a
    clear i
    
    %Save ROI tseries for next use
    if strcmp(o.YesorNo,'n')
        
        myROIi = [];
        fullroi = [o.myROIpath o.myROIname{ROIi} '.mat']; %roi with path
        myROIi = loadROITSeries(v,fullroi);        %load (very slow)
        mkdir(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])    %save once to reload
        cd(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])
        save(['tSeries' o.myGroup o.myROIname{ROIi} date],'myROIi')
        cd ..
        
    end
    
    %otherwise load roi tseries
    if strcmp(o.YesorNo,'y')
        cd(o.myPath)
        thistSeries = dir(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/*tSeries' o.myGroup o.myROIname{ROIi} '*']);
        cd(o.myPath)
        myROIi = importdata(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/' thistSeries.name]);
    end
        
    %start analyses
    fprintf('%s %s %s \n','---------- ', o.myROIname{ROIi},' -------------')
    c.myROIinfo = myROIi;
    
    %----------------------------------------------------
    %Bar accuracy for a specified period after event onset
    %-----------------------------------------------------
    if strcmp(o.myAnalysis,'accAtTime')
                
        fprintf('%s \n','------- Now decoding from specified period -------')
        c.blockLen = o.volsToClass(ROIi,2) - o.volsToClass(ROIi,1) + 1;
        c.startLag = o.volsToClass(ROIi,1);        
        c.startLag = ['startLag=' num2str(c.startLag)]; %get startlag
        c.blockLen = ['blockLen=' num2str(c.blockLen)]; %get block length 
        c.myROIinfo = getInstances(v,c.myROIinfo,o.stimvol,'n=inf',c.startLag,c.blockLen); %Create instances by averaging vols in a period after stim. onse
        c.myROIinfo = c.myROIinfo{1};
        c.i         = c.myROIinfo.classify.instances;           %Bold average over vols
        c.myClasf.raw.fullSg = leaveOneOut(c.i);                %RAW data
        [c.izsc,c.pStgs] = preprocessInstances(c.i,'zscore=1'); %z-scored data
        c.myClasf.Zsc.fullSg = leaveOneOut(c.izsc);
        
        %Map of voxels weight in classification
        %--------------------------------------
        %o.weightMap = getClassifierWeightMap(v,o.myROI,o.stimvol);
        %mrDispOverlay(o.weightMap,o.myROI.scanNum,o.myROI.o.groupNum,v,...
        %     'colormapType','normal','cmap',splitcmap,...
        %     'range=[-9 9]','clip=[-9 9]');
        
        %bar
        figure('color','w')
        
        %accuracy from raw
        bar(1,c.myClasf.raw.fullSg.correct,'facecolor',[.4 .4 .4])
        hold on;SLerrorbar(1, c.myClasf.raw.fullSg.correct, 'yError',c.myClasf.raw.fullSg.correctSTE,...
            ['Color=[',num2str([0 0 0]),']'],'MarkerSize=1','linewidth',2,'linesmoothing','on');
        
        bar(2,c.myClasf.Zsc.fullSg.correct,'facecolor',[.4 .4 .4])
        hold on;SLerrorbar(2, c.myClasf.Zsc.fullSg.correct, 'yError',c.myClasf.Zsc.fullSg.correctSTE,...
            ['Color=[',num2str([0 0 0]),']'],'MarkerSize=1','linewidth',2,'linesmoothing','on');
        
        %chance
        o.nClasses = length(o.stimvol);
        hline(1/o.nClasses)
        xlim([0 3])
        set(gca,'xtick',[1 2],'xticklabel',{'Raw','zscored'},'fontsize',12)
        xlabel('BOLD response preprocessing')
        ylabel('Decoding accuracy (% correct)')
        nVx = size(c.i{1},2);
        for i = 1 : o.nClasses
            o.repByClass(i) = size(o.stimvol{i},2);
        end
        title([o.myROIname{ROIi} ' - ' o.myClass ' (' num2str(o.nClasses) ' classes) - ' num2str(nVx) ' voxels - ' num2str(min(o.repByClass)) ' reps min per class'])
        box off
    end
    
    %------------------------------------------------------
    %Plot accuracy for contiguous pairs of vols after event
    %------------------------------------------------------
    %Create instances by averaging volumes in a defined period after stimulus
    %onset. Set lag in volumes and length of block data  in vol. over which to
    %average data into instances.
    %cover a up to 20 sec after event onset with time steps of 2 vols.
    %check analysis
    if strcmp(o.myAnalysis,'accByTime')
        
        fprintf('%s \n','------------------- Now decoding at contiguous pairs of vols -----------------')
        c.tr = 0.5;
        for i = 1 : round(20/c.tr)
            stLags{i} = ['startLag=' num2str(i)];
        end
        blockLens   = 'blockLen=2';
        myROI       = c.myROIinfo;
        stimvol     = o.stimvol;
        lengthsteps = numel(stLags);
        
        for j = 1 : numel(stLags)
            
            %classify raw instances
            %----------------------
            fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Decoding from raw Bold....')
            myROIj{j}     = getInstances(v,myROI,stimvol,'n=inf',stLags{j},blockLens);
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
            fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Decoding from z-scored Bold....')
            [izscj{j},pStgsj{j}] = preprocessInstances(ij{j},'zscore=1');
            ZscPerperiods{j}     = leaveOneOut(izscj{j});
            
            %print
            fprintf('%12s %i %s \n','(slfMRIclassifAnalSensoryAreas) ',round(j/lengthsteps *100),'% complete')
            myROIj = [];
        end
        
        %output
        c.stLags    = stLags;
        c.blockLens = blockLens;
        c.ij        = ij;
        
        %---------------------------
        %Plot accuracy by time steps
        %---------------------------
        figure('color','w')
        
        %get accuracies
        for i = 1:lengthsteps
            c.myClasf.raw.perperiods.Accy(i)   = perperiods{i}.correct;
            c.myClasf.raw.perperiods.AccySTE(i)= perperiods{i}.correctSTE;
            c.myClasf.Zsc.perperiods.Accy(i)   = ZscPerperiods{i}.correct;
            c.myClasf.Zsc.perperiods.AccySTE(i)= ZscPerperiods{i}.correctSTE;
        end
        
        %plot accuracy from raw instances
        SLerrorbar([1:lengthsteps]*c.tr, c.myClasf.raw.perperiods.Accy, 'yError',c.myClasf.raw.perperiods.AccySTE,...
            ['Color=[',num2str([.8 .8 .8]),']'],'MarkerSize=7','linewidth',2,'linesmoothing','on');
        
        %plot accuracy from zscored instances
        hold on; SLerrorbar([1:lengthsteps]*c.tr, c.myClasf.Zsc.perperiods.Accy, 'yError',c.myClasf.Zsc.perperiods.AccySTE,...
            ['Color=[',num2str([0 0 0]),']'],'MarkerSize=7','linewidth',2,'linesmoothing','on');
        
        legend('raw Bold (grey)','z-scored Bold (black)')
        
        %full backup
        c.perperiods{i}.correct = perperiods{i}.correct;
        c.perperiods{i}.correctSTE = perperiods{i}.correctSTE;
        c.ZscPerperiods{i}.correct = ZscPerperiods{i}.correct;
        c.ZscPerperiods{i}.correctSTE = ZscPerperiods{i}.correctSTE;
        
        %get task segments and add them to plot
        a = viewGet(v,'stimfile');
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
        xlim([min(TrialEventTimes)-1 max([TrialEventTimes lengthsteps*c.tr])])
        for i = 1 : length(seglens)+1
            vline(TrialEventTimes(i))
        end
        
        % with info
        o.nClasses = length(stimvol);
        randAccy = 1./o.nClasses;
        hline(randAccy)
        nVx = size(ij{1},2);
        for i = 1 : o.nClasses
            o.repByClass(i) = size(stimvol{i},2);
        end
        [~,c.ROIname] = fileparts(o.myROIname{ROIi});
        title([c.ROIname{ROIi} ' - ' o.myClass ' (' num2str(o.nClasses) ' classes) - ' num2str(nVx) ' voxels - ' num2str(min(o.repByClass)) ' reps min per class'])
        ylim([randAccy+[ -0.2 .2]])
        xlabel('Time after event (sec)')
        ylabel('Accy (% correct)')
        box off
        
        %-------
        %Summary
        %-------
        fprintf('-------- Results classif. (full period) -------------- \n')
        fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Accuracy Raw data     : ',c.myClasf.raw.fullSg.correct)
        fprintf('%s %i \n \n','(slfMRIclassifAnalSensoryAreas) Accuracy zscored data : ',c.myClasf.Zsc.fullSg.correct)
        
        fprintf('-------- Results classif (per period) --------------- \n')
        fprintf('%s %.2f \n','(slfMRIclassifAnalSensoryAreas) Accuracy Raw data     : ',c.myClasf.raw.perperiods.Accy)
        fprintf('%s %.2f \n \n','(slfMRIclassifAnalSensoryAreas) Accuracy zscored data : ',c.myClasf.Zsc.perperiods.Accy)
    end
        
    %backup
    cd(o.myPath)
    mkdir('slAnalyses'); cd slAnalyses
    mkdir(o.myGroup) ; cd(o.myGroup)
    mkdir('classif') ; cd classif
    mkdir(o.myClass) ; cd(o.myClass)
    mkdir(o.myAnalysis) ; cd(o.myAnalysis)
    mkdir(o.myROIname{ROIi}) ; cd(o.myROIname{ROIi}) 
    save(['Classif' o.myClass o.myROIname{ROIi} date '.mat'],'o','c')
    saveas(gcf, ['Classif' o.myClass o.myROIname{ROIi} '.fig'], 'fig')
    o.duration = toc;
   
end
