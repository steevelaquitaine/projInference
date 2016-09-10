
%slfMRIclassifAnalSensoryAreas2.m
%
%    Author: steeve laquitaine
%      date: 150922
%   purpose: decoding motion direction or coherence from ROI voxels Bold with Fisher linear discriminant
%
%     usage:
%                o.myPath = '~/dataold/datafMRI/sltaskdotdirfmri05/s02520150814'
%               o.myGroup = 'Concatenation';
%        o.myGroupScanNum = 1;
%              o.scanList = 1:18;                       %for concatenation if needed
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
%
%Analyses:
%                                            'accuracyByTime': plot accuracy by time averaged over contiguous pairs of vols after event
%           'accuracyAtTime',[1 2; 1 2; 1 2; 1 2; 1 2; 1 2]) : plot accuracy for instances averaged over vols 1 and 2 after event for many rois.
%                                                              [Start end;] are input for each Roi analysed.
%                                         'decodeF1forEachF2':
%
%
%
%
%
%note:
%
%   Duration : about 36 min for 18 x 614 vols.
%
%
%Organize output folders and data....currently a mess !!!

%
function [o,c,v] = slfMRIclassifAnalSensoryAreas2(o,varargin)

%-------
%ANALYSES
%--------
o.myAnalysis = [];

if sum(strcmp(varargin,'accuracyByTime'))
    %accuracy time course
    o.myAnalysis = {'accByTime'};
elseif sum(strcmp(varargin,'accuracyAtTime'))
    %accuracy at specific time
    o.myAnalysis = {'accAtTime'};
    pos = find(strcmp(varargin,'accuracyAtTime')) + 1;
    o.volsToClass = varargin{pos};                        %[start end] vols to average for classification
end
if sum(strcmp(varargin,'r2cutoff'))
    %accuracy at specific time
    o.myAnalysis{end+1} = 'r2cutoff';
    pos = find(strcmp(varargin,'r2cutoff')) + 1;
    o.r2thres = varargin{pos};
end
%leave-one instance out
if sum(strcmp(varargin,'leaveOneOut'))
    o.myAnalysis{end+1} = 'leaveOneOut';
end
%leave-one k mean instance out
%group instances into K average instances by class and leave-one-out
%classification. The idea is to average out past trial Bold confound.
if sum(strcmp(varargin,'slleaveOneOfkmeanInstanceOut'))
    o.myAnalysis{end+1} = 'slleaveOneOfkmeanInstanceOut';
    pos = find(strcmp(varargin,'slleaveOneOfkmeanInstanceOut')) + 1;
    o.kInstances = varargin{pos};
end
%regress out a variable
if any(strcmp(varargin,'regOutVar2'))
    o.myAnalysis{end+1} = 'regOutVar2';
    pos = find(strcmp(varargin,'regOutVar2')) + 1;
    o.regOutVar2 = varargin{pos};
end

%launch clean mrLoadRet
tic
cd(o.myPath)
v = mrLoadRet([]);

%--------
%Sum info
%--------
fprintf('--------------------- Data info ------------------------------ \n')
nRois = length(o.myROIname);
for i = 1 : nRois
    fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) ROI: ',o.myROIname{i})
end
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task: ',o.taskNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task segment: ',o.segmentNum)

fprintf('--------------------- Classification info -------------------- \n')
fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Fisher linear discriminant')
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Class:',o.myClass)
fprintf('-------------------------------------------------------------- \n')

%group & scan
v = viewSet(getMLRView,'curGroup',o.myGroup);
v = viewSet(getMLRView,'curScan',o.myScan);
refreshMLRDisplay(v); %update display

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

%if variable is regressed out get its stimvol
if any(strcmp(o.myAnalysis,'regOutVar2'))
  fprintf('%s %s %s %s \n','(slfMRIclassifAnalSensoryAreas2) Variable',o.regOutVar2,'will be regress-out before classification of variable',o.myClass)
  fprintf('%s %s %s %s \n','(slfMRIclassifAnalSensoryAreas2) Getting its stimvols')
  [o.stimvol2,o.taskCond2] = getStimvol(v,o.regOutVar2,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
end 

%----------------------------------------------
% Run / Load _All_ er & choose vox with high r2
%----------------------------------------------
%Check if er. already exists
if ~isempty(dir([o.myGroup '/' 'erAnal*']))
    fprintf('\n %s \n \n','-------------------- README -------------------')
    YorN = [];
    while ~strcmp(YorN,'y') && ~strcmp(YorN,'n')
      YorN = input(['!! README !! er. takes a while for long scans (e.g., concatenation)',...
		    '. If you"ve run er LOAD IT !! Otherwise make sure that you set',...
		    'maxBlocksize in Edit - preferences to a reasonable low value',...
		    '(e.g., 30GB) Do you want to load an already existing "erAnal" ',...
		    'for speed ? [y/n] :'],'s');      
    end
end

%if er. exists load
if strcmp(YorN,'y')
    fprintf('\n %s \n \n','------------ Loading existing erAnal ----------')
    v = loadAnalysis(v);
elseif strcmp(YorN,'n')
    %run er
    fprintf('\n %s \n \n',' ------------ Running er analysis ------------ ')
    
    %init default params
    [~,params] = eventRelated(v,[],'justGetParams=1','defaultParams=1',...
        'scanList',o.myScan);
    fprintf('%s %s %s \n','(slfMRIerAnal)','Available variables for er analysis are :',params.scanParams{1}.paramInfo{7}{3})
    params.applyFiltering                                          = 1;           %filtering
    params.scanParams{o.myScan}.stimvolVarInfo.varname     = '_all_';     %er variable
    fprintf('%s \n','(slfMRIerAnal) Running er for _all_')
    params.scanParams{o.myScan}.stimvolVarInfo.taskNum     = o.taskNum;   %task
    params.scanParams{o.myScan}.stimvolVarInfo.segmentNum  = o.segmentNum;
    %params.scanParams{o.myScan}.varname                    = o.myVar;     %er variable
    %params.scanParams{o.myScan}.taskNum                    = o.taskNum;   %task
    %params.scanParams{o.myScan}.segmentNum                 = o.segmentNum;
    [v,d] = eventRelated(v,params); %r2 map
else
  keyboard
end

%get r2
o.r2All      = viewGet(v,'overlayData',o.myScan);
o.scanDims   = viewGet(v,'scanDims'); %N by N by # slice

%------------------------- load ROITseries -------------------------------
%This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
%GET rid of it asap. Prob is er analysis doesn't load and store the
%entire volume but load slice by slice. Thus cannot be used here to get
%the ROI TSeries. The fastest way so far is to save the ROI TSeries once
%run once and reload them instead of loadROITSeries them.
%e.myROI.scanCoords = getROICoordinates(v,myROIname);
o.YesorNo = [];
while ~strcmp(o.YesorNo,'y') && ~strcmp(o.YesorNo,'n')
   o.YesorNo = input(['!! README !! Loading ROI tseries takes a while for long'...
		      ' scans (e.g., concatenation). Do you want to load existing ROItseries ? [y/n]'],'s');
end
 
%if no, load roi tseries
o.nRois = length(o.myROIname);

%-----------------
%CLASSIFY from ROI
%-----------------
%repeat analyses by rois
for ROIi = 1 : nRois
    
    %first cleanup
    keep = {'v','o','ROIi','keep','a','i'}; a = whos;
    for i = 1 : length(a)
        if ~sum(strcmp(a(i).name,keep))>0
            clear(a(i).name)
        end
    end
    clear keep
    clear a
    clear i
    
    %Save ROI or/and use existing ROI tseries
    if strcmp(o.YesorNo,'n')
        myROIi = [];
        fullroi = [o.myROIpath o.myROIname{ROIi} '.mat']; %roi with path
        myROIi = loadROITSeries(v,fullroi);               %load (very slow)
        mkdir(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])    %save once to reload
        cd(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])
        save(['tSeries' o.myGroup o.myROIname{ROIi} date],'myROIi')
        cd ..
    elseif strcmp(o.YesorNo,'y')
        %otherwise load roi tseries
        cd(o.myPath)
        thistSeries = dir(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/*tSeries' o.myGroup o.myROIname{ROIi} '*']);
        cd(o.myPath)
        myROIi = importdata(['slAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/' thistSeries.name]);
    end
         
    %------------- simulate ROITSeries to check code (to comment) ---------
    %fprintf('%s %s %s \n','(slfMRIerAnal)','Simulating ROI TSeries to check code')
    %myROIi = slSimROITSeries(myROIi,o);
    %------------------------------------------------------------------------
    
    %start analyses
    fprintf('\n %s %s %s \n','---------- ', o.myROIname{ROIi},' -------------')
    c.myROIinfo = myROIi;
    
    %r2 cutoff
    %---------
    if any(strcmp(o.myAnalysis,'r2cutoff'))
        r2cutoff = ['r2cutoff=' num2str(o.r2thres)];
        c.myROIinfo = slSortROIvoxels(myROIi,r2cutoff);
    end
    
    %----------------------------------------------------
    %Bar accuracy for a specified period after event onset
    %-----------------------------------------------------
    if any(strcmp(o.myAnalysis,'accAtTime'))
        
        %get instances (Bold responses) for variable 1
        fprintf('\n %s \n','------- Now decoding from specified time period -------')
        c.blockLen = o.volsToClass(ROIi,2) - o.volsToClass(ROIi,1) + 1;
        c.startLag = o.volsToClass(ROIi,1);
        c.startLag = ['startLag=' num2str(c.startLag)];
        c.blockLen = ['blockLen=' num2str(c.blockLen)];		
        c.myROIinfo1 = getInstances(v,c.myROIinfo,o.stimvol,'n=inf', c.startLag,c.blockLen); %average blocklen vols after stim. onset
        c.myROIinfo1 = c.myROIinfo1{1};
        c.i          = c.myROIinfo1.classify.instances;      %Bold average over vols
        c.stimvol1   = c.myROIinfo1.classify.instanceVol;

        %regress variable 2 weights out of instances
	%the weigths of variable 1 classes and variable 2 classes are evaluated with linear regression
	%Instances are then recalculated by summing the weights of variable 1, the regression constant and 
	%error term but excluding variable 2 classes weights.
        c.info2    = getInstances(v,c.myROIinfo,o.stimvol2,'n=inf',c.startLag,c.blockLen); 	
	c.i2       = c.info2{1}.classify.instances;
	c.stimvol2 = c.info2{1}.classify.instanceVol;
        if any(strcmp(o.myAnalysis,'regOutVar2'))
	  [c.i,o.taskCond,o.stimVol] = slregOutVarEffectOnInstance(c.i,o.taskCond,c.stimvol1,c.i2,o.taskCond2,c.stimvol2);
	end
	
        %leave-one-instance out fisher classification
        if any(strcmp(o.myAnalysis,'leaveOneOut'))
	    fprintf('%s \n','(leaveOneOut) Running leaveOneOut.')
	    %raw
            c.myClasf.raw.fullSg = leaveOneOut(c.i);               
	    %zscore
            [c.izsc,c.pStgs] = preprocessInstances(c.i,'zscore=1');
            c.myClasf.Zsc.fullSg = leaveOneOut(c.izsc);
        end
	
        %leave-one Kmean instance out fisher classification
        if any(strcmp(o.myAnalysis,'slleaveOneOfkmeanInstanceOut'))
            fprintf('%s \n','(slleaveOneOfkmeanInstanceOut) Running slleaveOneOfkmeanInstanceOut.')
	    %raw
            c.myClasf.raw.fullSg = slleaveOneOfkmeanInstanceOut(c.i,o.kInstances);
	    %szcore
            [c.izsc,c.pStgs] = preprocessInstances(c.i,'zscore=1'); 
            c.myClasf.Zsc.fullSg = slleaveOneOfkmeanInstanceOut(c.izsc,o.kInstances);
        end
        
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
    if any(strcmp(o.myAnalysis,'accByTime'))
        
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
            
            %leave-one-out instance out
            if any(strcmp(o.myAnalysis,'leaveOneOut'))
                perperiods{j} = leaveOneOut(ij{j});                %RAW data
                [izscj{j},pStgsj{j}] = preprocessInstances(ij{j},'zscore=1'); %z-scored data
                ZscPerperiods{j} = leaveOneOut(c.izsc);
            end
            %leave-one Kmean instance out
            if any(strcmp(o.myAnalysis,'slleaveOneOfkmeanInstanceOut'))
                fprintf('%s \n','(slleaveOneOfkmeanInstanceOut) Running slleaveOneOfkmeanInstanceOut.')
                perperiods{j} = slleaveOneOfkmeanInstanceOut(ij{j},o.kInstances);
                [izscj{j},pStgsj{j}] = preprocessInstances(ij{j},'zscore=1');
                ZscPerperiods{j}  = slleaveOneOfkmeanInstanceOut(c.izsc,o.kInstances);
            end
                        
            %plot MAP of voxels weight
            %o.weightMap = getClassifierWeightMap(v,o.myROI,stimvol);
            
            %remove current analysis and look at the classifier's weightMap
            % % - view menu - remove Analysis
            %mrDispOverlay(weightMap,myROI.scanNum,myROI.groupNum,getMLRView,...
            %     'colormapType','normal','cmap',splitcmap,...
            %     'range=[-9 9]','clip=[-9 9]');
            
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
    mkdir([o.myAnalysis{:}]) ; cd([o.myAnalysis{:}])
    mkdir(o.myROIname{ROIi}) ; cd(o.myROIname{ROIi})
    save(['Classif' o.myClass o.myROIname{ROIi} date '.mat'],'o','c')
    saveas(gcf, ['Classif' o.myClass o.myROIname{ROIi} '.fig'], 'fig')
    o.myAnalPath = pwd;
    o.duration = toc;
    
end



