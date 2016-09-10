
%slfmriClassSensStckSess2.m
%
%    Author: steeve laquitaine
%      date: 150922
%   purpose: decoding motion direction or coherence from ROI voxels Bold with Fisher linear discriminant
%
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
%
%to do : # of instances are not consistent per session. how to handle that
%? !!!!!!!!
%
function [si,o,c] = slfmriClassSensStckSess2(o,myClass,varargin)

o.myClass = myClass;

%-------
%ANALYSES
%--------
mrQuit
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

if slIsInput(varargin,'CalcInstances')
    %------------------------- Get aligned instances for classif ----------------------
    %loop over sessions
    tic
    nSess = length(o.sessPath);
    
    for si = 1 : nSess
        
        %go to data
        cd(o.sessPath{si})
        h = mrLoadRet([]);
        v{si} = h;
        
        fprintf('%s \n', ['Session:' o.sessPath{si} ' ---------------------------'])
        
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
        v{si} = viewSet(getMLRView,'curGroup',o.myGroup);
        v{si} = viewSet(getMLRView,'curScan',o.myScan);
        refreshMLRDisplay(v{si}); %update display
        
        %----------
        %INFO CHECK
        %----------
        %Vols to decode class at event (Nclasses cells of 1 by Ntrialrep vol #)
        [d{si}, o.taskCond] = getStimvol(v{si},o.myClass,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
        if isfield(d{si},'stimvol')
            o.stimvol = d{si}.stimvol;
        else
            o.stimvol = d{si};
        end
        fprintf('%s %i \n','---------------- Please check that those info are correct -----------------------------')
        fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Current group: ', viewGet(v{si},'curGroup'))
        fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Current scan: ', viewGet(v{si},'curScan'))
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
            [o.stimvol2,o.taskCond2] = getStimvol(v{si},o.regOutVar2,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
        end
        
        %------------------------- load ROITseries -------------------------------
        %This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
        %GET rid of it asap. Prob is er analysis doesn't load and store the
        %entire volume but load slice by slice. Thus cannot be used here to get
        %the ROI TSeries. The fastest way so far is to save the ROI TSeries once
        %run once and reload them instead of loadROITSeries them.
        %e.myROI.scanCoords = getROICoordinates(v,myROIname);
        o.YesorNo = [];
        o.YesorNo = 'y';
        while ~strcmp(o.YesorNo,'y') && ~strcmp(o.YesorNo,'n')
            fprintf('%s \n', '------------------- Getting ROI Tseries ------------------')
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
            t1=tic;
            %first cleanup
            keep = {'s','si','v','o','ROIi','keep','a','i'}; a = whos;
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
                myROIi = loadROITSeries(v{si},fullroi);               %load (very slow)
                mkdir([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])    %save once to reload
                cd([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])
                save(['tSeries' o.myGroup o.myROIname{ROIi} date],'myROIi')
                cd ..
            elseif strcmp(o.YesorNo,'y')
                %otherwise load roi tseries
                cd(o.sessPath{si})
                thistSeries = dir([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/*tSeries' o.myGroup o.myROIname{ROIi} '*']);
                cd(o.sessPath{si})
                myROIi = importdata([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/' thistSeries.name]);
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
                fprintf('\n %s \n','------- Now averaging instances from specified time period -------')
                c.blockLen = o.volsToClass(ROIi,2) - o.volsToClass(ROIi,1) + 1;
                c.startLag = o.volsToClass(ROIi,1);
                c.startLag = ['startLag=' num2str(c.startLag)];
                c.blockLen = ['blockLen=' num2str(c.blockLen)];
                c.myROIinfo1 = getInstances(v{si},c.myROIinfo,o.stimvol,'n=inf', c.startLag,c.blockLen); %average blocklen vols after stim. onset
                keyboard
                c.myROIinfo1 = c.myROIinfo1{1};
                c.i          = c.myROIinfo1.classify.instances;      %Bold average over vols
                c.stimvol1   = c.myROIinfo1.classify.instanceVol;
                
                %store session
                %raw
                s.i{ROIi}{si} = c.i;
                %zscored
                [c.izsc,c.pStgs] = preprocessInstances(c.i,'zscore=1');
                s.izsc{ROIi}{si} = c.izsc;
                
            end
        end
        mrQuit
    end
    
    %backup data
    cd([o.stckPath 'slStckAnalyses/' o.myGroup])
    if ~isdir('classif'); mkdir('classif');end
    cd('classif')
    if ~isdir(o.myClass); mkdir( o.myClass);end
    cd(o.myClass)
    if ~isdir([o.myAnalysis{:}]); mkdir([o.myAnalysis{:}]);end
    cd([o.myAnalysis{:}])
    save(['ClassifStckSessData' o.myClass '_allRois_' date '.mat'],'o','c','s')
    
    
    %------------------------------- Classification ------------------------
    %load existing instance data
elseif slIsInput(varargin,'loadSavedInstances')
    cd([o.stckPath 'slStckAnalyses/' o.myGroup])
    cd('classif')
    cd(o.myClass)
    cd([o.myAnalysis{:}])
    l = dir(['ClassifStckSessData' o.myClass '_allRois_*']);
    if isempty(l)
        slPrintfStr('slfmriClassSensStckSess',[' ClassifStckSessData' o.myClass '_allRois_..not found...'])
        keyboard
    else
        load(l.name)
    end
else
    slPrintfStr('slfmriClassSensStckSess','Please input either ')
    slPrintfStr('slfmriClassSensStckSess','- loadSavedInstances : to load existing instances')
    slPrintfStr('slfmriClassSensStckSess','- CalcInstances : to calculate instances')
    keyboard
end

%classify
nVar = length(o.taskCond);
t1 = tic;

%loop over roi
for roii = 1 : o.nRois
    
    fprintf('%s \n','------------------- Classification ------------------')
    
    %path
    cd([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.myClass '/' [o.myAnalysis{:}]])
    
    %init
    c.myClasf.raw.fullSg = [];
    c.myClasf.Zsc.fullSg = [];
    
    %loop over variables
    for vi = 1 : nVar
        %loop over sessions
        %stki{roii}{vi} = [];
        Aligni{roii}{vi} = [];
        Aligniz{roii}{vi} = [];
        for si = 1 :  length(o.sessPath)
            if vi > length(s.i{roii}{si})
                s.i{roii}{si}{vi}={};
            end
            if vi > length(s.izsc{roii}{si})
                s.izsc{roii}{si}{vi}={};
                isz=[];
            end
            is = s.i{roii}{si}{vi};
            isz = s.izsc{roii}{si}{vi};
            if isempty(s.i{roii}{si}{vi})
                is =[];
            end
            if isempty(s.izsc{roii}{si}{vi})
                isz =[];
            end
            %Instances (col) by vox (rows) stacked across sessions
            %sorted by classes (cells)
            Aligni{roii}{vi} = [Aligni{roii}{vi} ; is]; %raw
            Aligniz{roii}{vi} = [Aligniz{roii}{vi} ; isz];%z-scored
        end
    end
    
    fprintf('%s \n',['(leaveOneOut) Running leaveOneOut for ' o.myROIname{roii}] )
    
    %------------ classification ----------------------
    %raw
    c.myClasf.raw.fullSg = leaveOneOut(Aligni{roii});
    %zscore
    c.myClasf.Zsc.fullSg = leaveOneOut(Aligniz{roii});
    
    %------------ classification (shuffled labels) ----------------------
    fprintf('%s \n','------------------- Classification (Shuffled) ------------------')
    
    %shuffle between classes
    %concat all instances
    AlgAll{roii} = cell2mat(Aligni{roii}');   %raw
    AlgAllz{roii} = cell2mat(Aligniz{roii}'); %z-sc
    %get their classes
    for vi = 1 : nVar
        ni(vi) = size(Aligni{roii}{vi},1);
    end
    niend = cumsum(ni);
    nist = [1 ni(1:end-1)+1];
    %initialize class-shuffling
    nSf = 10; %default: nSf=50: about 3h (10x is 16 min)
    mA = repmat(AlgAll(roii),nSf,1);
    mAz = repmat(AlgAllz(roii),nSf,1);
    nistAll = repmat(nist,nSf,1);
    niendAll = repmat(niend,nSf,1);
    %create n shuffled-class data
    for ip = 1 : nSf
        %shuffle class associated with instance
        shf = randperm(size(mA{ip},1));
        AlgAllsh = mA{ip}(shf,:);   %raw : Ntotalinst by Nvox
        AlgAllshz = mAz{ip}(shf,:); %z-scored
        for vi = 1 : nVar
            AlgSh{ip}{vi} = AlgAllsh(nistAll(ip,vi):niendAll(ip,vi),:);
            AlgShz{ip}{vi} = AlgAllshz(nistAll(ip,vi):niendAll(ip,vi),:);
        end
    end
    %classification of the n shuffled dataset
    parfor ip = 1 : nSf
        cip(ip) = leaveOneOut(AlgSh{ip});  %raw
        czip(ip) = leaveOneOut(AlgShz{ip});%zscore
        fprintf('%s ','---------------------------')
    end
    c.myClasf.raw.fullSgShuf = cip;
    c.myClasf.Zsc.fullSgShuf = czip;
    
    fprintf('%s \n','------------------- Classification (balanced) ------------------')
    %Check that classification training is not affected by
    %the dataset imbalance. remove instances from the
    %class with most instances to get equal instance proportions
    %between all classes
    cb = classifBal(Aligni{roii},Aligniz{roii});
    c.myClasf.raw.fullSgShufBal = cb.myClasf.raw.fullSgShufBal;
    c.myClasf.Zsc.fullSgShufBal = cb.myClasf.Zsc.fullSgShufBal;
        
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    o.myAnalPath2{roii} = pwd;
    o.duration2 = toc(t1);
    save(['ClassifStckSess' o.myClass '_' o.myROIname{roii} '_' date '.mat'],'o','c','s')
    
    %cleanup
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShuf');%shuffled
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShuf');%shuffled
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShufBal');%shuffled
    c.myClasf.raw = rmfield(c.myClasf.Zsc,'fullSgShufBal');%shuffled       
end

%Balance data set then classify
%input: cells of instances per class (raw and zscored)
function cb = classifBal(rawI,zI)

%------------------ balance classes -----------------------
%get #of instance by class
nClass = length(rawI);
for i = 1 : nClass
    nI(i) = size(rawI{i},1);
end

%get minimum # of instance
minnI = min(nI);

%set all classes # of
%instances to the minimum
for i = 1 : nClass
    rawIbal{i} = rawI{i}(1:minnI,:); %raw
    zIbal{i}   = zI{i}(1:minnI,:);   %zscored
end

%------------------ shuffle and classify -----------------------
%shuffle classes
%concat all instances
AlgAll = cell2mat(rawIbal');%raw
AlgAllz = cell2mat(zIbal'); %z-sc
%mark instance classes
ni = repmat(minnI,1,nClass);
niend = cumsum(ni);
nist = [1 ni(1:end-1)+1];
%init class-shuffling
nSf = 10; %default: nSf=50: about 3h (10x is 16 min)
mA = repmat({AlgAll},nSf,1);
mAz = repmat({AlgAllz},nSf,1);
nistAll = repmat(nist,nSf,1);
niendAll = repmat(niend,nSf,1);
%create n shuffled-class data
for ip = 1 : nSf
    %shuffle class associated with instance
    shf = randperm(size(mA{ip},1));
    AlgAllsh = mA{ip}(shf,:);   %raw : Ntotalinst by Nvox
    AlgAllshz = mAz{ip}(shf,:); %z-scored
    for vi = 1 : nClass
        AlgSh{ip}{vi} = AlgAllsh(nistAll(ip,vi):niendAll(ip,vi),:);
        AlgShz{ip}{vi} = AlgAllshz(nistAll(ip,vi):niendAll(ip,vi),:);
    end
end
%classify each shuffled dataset
parfor ip = 1 : nSf
    cip(ip) = leaveOneOut(AlgSh{ip});  %raw
    czip(ip) = leaveOneOut(AlgShz{ip});%zscore
    fprintf('%s ','---------------------------')
end
cb.myClasf.raw.fullSgShufBal = cip;
cb.myClasf.Zsc.fullSgShufBal = czip;









