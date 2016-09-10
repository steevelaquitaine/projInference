
%slfmriClassification.m
%
%    Author: steeve laquitaine
%      date: 151217
%   purpose: decoding motion direction or coherence from ROI voxels Bold with Fisher linear discriminant
%
%inputs:
%
%       'CalcInstances'(or 'loadSavedInstances')
%       'slReSampInst': re-sample instances when few only
%                       set random values for missing classes. Do not run
%                       when classes are missing !!
%          'nPerm=100': to get chance level form 100 dataset of permuted 
%                       classes

function slfmriClassification(o,myClass,varargin)

%duration
t0=tic;

%getArgs
nPerm = [];
getArgs(varargin,{'nPerm'});

%get classes
o.myClass = myClass;

%get analysis
o.myAnalysis = [];
o = slgetAnalysis(o,varargin);

%get session-stacked instances
fprintf('%s \n','---- Get session-stacked instances ---')

%get instances
%calculate
if slIsInput(varargin,'CalcInstances')
    %get instances (roi,session)
    [o,si,s,c] = slCalcStackedInstances(o);
    %backup instances in structured directory
    [o,c,s] = slSavedInStrDirectory(o,c,s);
    
    %or load
elseif slIsInput(varargin,'loadSavedInstances')
    %go to path
    %convert cell class to string
    o.savedClass = o.myClass;
    while iscell(o.savedClass)
        o.savedClass = [o.savedClass{:};];
    end
    cd([o.stckPath 'slStckAnalyses/' o.myGroup 'classif' o.savedClass [o.myAnalysis{:}]])
    %check file exists
    l = dir(['ClassifStckSessData' o.myClass '_allRois_*']);
    if isempty(l)
        slPrintfStr('slfmriClassSensStckSess',...
            [' ClassifStckSessData' o.savedClass '_allRois_..not found...'])
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

%classify by roi
nVar = length(o.taskCond);
for roii = 1 : o.nRois
    fprintf('%s \n',['(leaveOneOut)' o.myROIname{roii}] )
    
    %path
    cd([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]])
    
    %init
    [Aligni,Aligniz,s] = slStackSessInst(o,s,nVar,roii);
    
    %!!!!TEST WARNING !!!!
    %----simulate instances for test
    %c=simInstances({'lV1'},nVar,[],repmat(25,nVar,1),'type=a2Dclusters');
    %Aligni=repmat({c{1}.classify.instances},o.nRois,1);
    %Aligniz=preprocessInstances(Aligni{1},'zscore=1');
    %Aligniz=repmat({Aligniz},o.nRois,1);
    %!!!!!!!!!!!!!! debugging ----------------
    
    %!!!!TEST WARNING !!!! resample instances to increase
    %size of dataset.just for test: might break
    %the training/test isolation requirement.
    %also I set to rand the classes when instances
    %that have no instances
    if slIsInput(varargin,'slReSampInst')
        [Aligni,Aligniz] = slReSampInst(Aligni,Aligniz);
    end
    %!!!!!!!!!!!!!! debugging ----------------
    
    keyboard
    
    
    
    
    
    fprintf('%s \n','------ Classification (actual data) --------')
    c.myClasf.raw.fullSg = [];
    c.myClasf.Zsc.fullSg = [];
    %case leaveoneout
    if slIsInput(varargin,'leaveOneOut')
        c.myClasf.raw.fullSg = leaveOneOut(Aligni{roii});%raw
        c.myClasf.Zsc.fullSg = leaveOneOut(Aligniz{roii});%zs
    end
    %case 10-fold
    if slIsInput(varargin,'kFoldCV')
        fprintf('%s \n','kFoldCV')
        fprintf('%s \n \n','-----------')
        c.myClasf.raw.fullSg = kFold(Aligni{roii},'numFolds=10');%raw
        c.myClasf.Zsc.fullSg = kFold(Aligniz{roii},'numFolds=10');%zs
    
    end
    
    fprintf('%s \n','---- (Classification at chance) permuted data -------')
    %classify permutated data ntimes
    %36 min (50 perms,1 ROI)
    %typically 30 is enough
    %(Markus Ojala, see Permutation test for stuyding classifer performance,
    %Journal of Machine learning, 2010)
    c.myClasf.raw.fullSgShuf=[];
    c.myClasf.Zsc.fullSgShuf=[];
    %case leaveOneOut
    if slIsInput(varargin,'leaveOneOut')        
        if ~isempty(nPerm)
            fprintf('%s %i %s \n','(slfmriClassification) Running ',nPerm,'permutations...')
            nPerm2 = ['nPerm=' num2str(nPerm)];
            c.myClasf.raw.fullSgShuf = leaveOneOutNpermut(Aligni{roii},'permutationUnBal=1',nPerm2);
            c.myClasf.Zsc.fullSgShuf = leaveOneOutNpermut(Aligni{roii},'permutationUnBal=1',nPerm2);
        end
    end
    %case 10-Fold
    if slIsInput(varargin,'kFoldCV')        
        if ~isempty(nPerm)
            fprintf('%s %i %s \n','(slfmriClassification) Running ',nPerm,'permutations...')
            nPerm2 = ['nPerm=' num2str(nPerm)];
            c.myClasf.raw.fullSgShuf = kFoldNpermut(Aligni{roii},nPerm2,'numFolds=10');
            c.myClasf.Zsc.fullSgShuf = kFoldNpermut(Aligni{roii},nPerm2,'numFolds=10');
        end
    end
    %fprintf('%s \n','------ Classification (balanced data) ---------')
    %uncomment. accuracy doesn't change much when we balance the dataset
    %balance data (removing) and classify
    c.myClasf.raw.fullSgShufBal = [];
    c.myClasf.Zsc.fullSgShufBal = [];
    %c.myClasf.raw.fullSgShufBal = leaveOneOut(Aligni{roii},'balancByRemovI=1');%raw
    %c.myClasf.Zsc.fullSgShufBal = leaveOneOut(Aligniz{roii},'balancByRemovI=1');%zscore
    
    %save results classification
    [o,c,s,roii] = slSaveClassiRes(o,c,s,roii);
    
    %duration
    duration=toc(t0);
    fprintf('%s %.2f \n','one more roi --------------',duration)
end

%back home
cd(o.rootpath)


%--------------- nested -------------
%get analysis
function o = slgetAnalysis(o,varargin)
%get arg
varargin = varargin{:};
%calculate accuracy each time step
if sum(strcmp(varargin,'accuracyByTime'))
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
%k-Fold
if sum(strcmp(varargin,'kFoldCV'))
    o.myAnalysis{end+1} = 'kFoldCV';
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
%store varargin
o.myvarg = varargin;

%display analysis info
function o = sldispAnalysisInfo(o)
fprintf('----- Analysis info -------- \n')
for i = 1 : o.nRois
    fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) ROI: ',o.myROIname{i})
end
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Base anat: ',o.anatFileName)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task: ',o.taskNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task phase: ',o.phaseNum)
fprintf('%s %i \n','(slfMRIclassifAnalSensoryAreas) Task segment: ',o.segmentNum)
fprintf('%s \n','(slfMRIclassifAnalSensoryAreas) Fisher linear discriminant')
%convert cell class to string
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
fprintf('%s %s \n','(slfMRIclassifAnalSensoryAreas) Class:',o.savedClass)



%calculate instances
function [o,si,s,c] = slCalcStackedInstances(o)

%mrLoadRet
tic

%each session
o.nSess = length(o.sessPath);
o.nRois = length(o.myROIname);
s = [];
for si = 1 : o.nSess
    
    %go to data
    cd(o.sessPath{si})
    
    %have to load each new session
    h=mrLoadRet([]);
    fprintf('%s \n', ['Session:' o.sessPath{si} ' ------'])
    
    %disp info
    o = sldispAnalysisInfo(o);
    
    %set group & scan
    v{si} = h;
    v{si} = viewSet(v{si},'curGroup',o.myGroup);
    v{si} = viewSet(v{si},'curScan',o.myScan);
    refreshMLRDisplay(v{si}); %update display
    
    %get each session stim vols
    o = slGetSessStimvol(o,v,si);
    
    %get each roi instances
    [o,v,si,s,c] = slgetROIinstances(o,v,si,s);
    mrQuit
end

%get session stimvol by class
function o = slGetSessStimvol(o,v,si)

%get stimvols by class
[d{si},o.taskCond] = getStimvol(v{si},o.myClass,'taskNum',o.taskNum,...
    'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);

%check if stimvol field exists
if isfield(d{si},'stimvol');
    o.stimvol = d{si}.stimvol;
else
    o.stimvol = d{si};
end

%session basic info
curG = viewGet(v{si},'curGroup');
curS = viewGet(v{si},'curScan');
maxvols = max([o.stimvol{:}]);
fprintf('%s %i \n','(slGetSessStimvol) Group: ',curG)
fprintf('%s %i \n','(slGetSessStimvol)  Scan: ',curS)
fprintf('%s %i \n','(slGetSessStimvol) Max # of vols: ',maxvols)
fprintf('%s \n','(slGetSessStimvol) Classes are: ')
for i = 1 : length(o.taskCond)
    disp(o.taskCond{i})
end

%if variable is regressed out get its stimvol
if any(strcmp(o.myAnalysis,'regOutVar2'))
    fprintf('%s %s %s %s \n','(slfMRIclassifAnalSensoryAreas2) Variable',o.regOutVar2,'will be regress-out before classification of variable',o.myClass)
    fprintf('%s %s %s %s \n','(slfMRIclassifAnalSensoryAreas2) Getting its stimvols')
    [o.stimvol2,o.taskCond2] = getStimvol(v{si},o.regOutVar2,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
end

%get session instances by roi
function [o,v,si,s,c] = slgetROIinstances(o,v,si,s)

%each roi
for ROIi = 1 : o.nRois
    
    %cleanup
    keep={'s','si','v','o','ROIi','keep','a','i'}; a=whos;
    for i=1:length(a)
        if ~sum(strcmp(a(i).name,keep))>0; clear(a(i).name); end
    end
    clear keep; clear a; clear i
    
    %re-loadROITseries (very long)
    if ~slIsInput(o.myvarg,'loadSavedROI')
        myROIi = slsaveROItseries(v,o,ROIi,si);
        %or load saved roi tseries
    elseif slIsInput(o.myvarg,'loadSavedROI')
        [thistSeries,myROIi] = slimportROItseries(o,si,ROIi);
    end
    
    %simulate ROITSeries----
    %fprintf('%s %s %s \n','(slfMRIerAnal)',...
    %'Simulating ROI TSeries to check code')
    %myROIi = slSimROITSeries(myROIi,o);
    %-----------------------------------
    
    %start analyses
    fprintf('\n %s %s %s \n','---- ', o.myROIname{ROIi},' ------')
    c.myROIinfo = myROIi;
    
    %r2 cutoff
    if any(strcmp(o.myAnalysis,'r2cutoff'))
        r2cutoff = ['r2cutoff=' num2str(o.r2thres)];
        c.myROIinfo = slSortROIvoxels(myROIi,r2cutoff);
    end
    
    %instances at a time
    if any(strcmp(o.myAnalysis,'accAtTime'))
        [s,v,o,c] = slgetInstancesAtTime(v,o,c,ROIi,si,s);
    end
end

%save ROI tseries in structured directory
function myROIi = slsaveROItseries(v,o,ROIi,si)
%roi with path
fullroi = [o.myROIpath o.myROIname{ROIi} '.mat']; 
%load ROITseries
%This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
%GET rid of it asap. Prob is er analysis doesn't load and store the
%entire volume but load slice by slice. Thus cannot be used here to get
%the ROI TSeries. The fastest way so far is to save the ROI TSeries once
%run once and reload them instead of loadROITSeries them.
%e.myROI.scanCoords = getROICoordinates(v,myROIname);
myROIi = loadROITSeries(v{si},fullroi); %load roi tseries (very slow)
%check roi
if isempty(myROIi)
    slPrintfStr('slfmriClassSensStckSess',' Could not find ROI. Change path ?')
    mrQuit; keyboard
end
%check dir exists
savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}];
if isdir(savedDir)
    mkdir([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])    %save once to reload
end
%backup
cd([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi}])
save(['tSeries' o.myGroup o.myROIname{ROIi} date],'myROIi')
cd ..

%load ROI tseries in workspace
function [thistSeries,myROIi] = slimportROItseries(o,si,ROIi)

cd(o.sessPath{si})
thistSeries = dir([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/*tSeries' o.myGroup o.myROIname{ROIi} '*']);
cd(o.sessPath{si})
myROIi = importdata([o.stckPath 'slStckAnalyses/' o.myGroup '/roiTseries/' o.myROIname{ROIi} '/' thistSeries.name]);

%get instances at time
function [s,v,o,c] = slgetInstancesAtTime(v,o,c,ROIi,si,s)

%get instances (Bold responses) for variable 1
fprintf('\n %s \n','------- Now averaging instances from specified time period -------')
%block length
c.blockLen = o.volsToClass(ROIi,2) - o.volsToClass(ROIi,1) + 1;
%lag to event
c.startLag = o.volsToClass(ROIi,1);
c.startLag = ['startLag=' num2str(c.startLag)];
c.blockLen = ['blockLen=' num2str(c.blockLen)];
%average blocklen vols after stim. onset
c.myROIinfo1 = getInstances(v{si},c.myROIinfo,o.stimvol,'n=inf',c.startLag,c.blockLen);
c.myROIinfo1 = c.myROIinfo1{1};
c.i          = c.myROIinfo1.classify.instances;      %Bold average over vols
c.stimvol1   = c.myROIinfo1.classify.instanceVol;

%store session
%raw
s.i{ROIi}{si} = c.i;
%zscored
[c.izsc,c.pStgs] = preprocessInstances(c.i,'zscore=1');
s.izsc{ROIi}{si} = c.izsc;

%stack instances across sessions
function [Aligni,Aligniz,s] = slStackSessInst(o,s,nVar,roii)

%loop over variables
for vi = 1 : nVar
    %collect instances each session
    Aligni{roii}{vi} = [];
    Aligniz{roii}{vi} = [];
    %stack instances over sessions
    for si = 1 :  length(o.sessPath)
        %deal with missing classes
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

%resample instances to increase
%sample size
function [Aligni,Aligniz] = slReSampInst(Aligni,Aligniz)

%get # vox
nClasses = length(Aligni{:});
nvox=[];
for cl = 1 : nClasses
    nvox = max([nvox size(Aligni{:}{cl})]);
end

nReSamp = 33;
for cl = 1 : nClasses
    nI = size(Aligni{:}{cl},1);
    %if no data
    if isempty(Aligni{:}{cl})
        Aligni{:}{cl} = rand(20,nvox);%raw
        Aligniz{:}{cl} = rand(20,nvox);%zscore
    else
        nReSampix = randsample(1:nI,nReSamp,'true');
        Aligni{:}{cl} = Aligni{:}{cl}(nReSampix,:);%raw
        Aligniz{:}{cl} = Aligniz{:}{cl}(nReSampix,:);%zscore
    end
end

%save stacked instances in structured directory
function [o,c,s] = slSavedInStrDirectory(o,c,s)
%structured directory
%convert cell class to string
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
%move to path
mkdir([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]]);
cd([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]]);
o.savedallRoisInstancesPath = pwd;
%set file name
filename = ['ClassifStckSessData' o.savedClass '_allRois_' date '.mat'];
%archive existing file
list = dir('*.mat');
if ~isempty(list)
    if isempty(dir('Archive')); mkdir Archive; end
    for i=1:length(list); movefile(list(i).name,'Archive/'); end
end
%save file
save(filename,'o','c','s')


%save classification data
function [o,c,s,roii] = slSaveClassiRes(o,c,s,roii)

%convert cell class to string
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
%backup
if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
cd(o.myROIname{roii})
o.myAnalPath2{roii} = pwd;
%archive existing file
list = dir('*.mat');
if ~isempty(list)
    if isempty(dir('Archive')); mkdir Archive; end
    for i=1:length(list); movefile(list(i).name,'Archive/'); end
end
save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '_' date '.mat'],'o','c','s')

%cleanup
if isfield(c.myClasf.raw,'fullSgShuf')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShuf');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShuf');
end
if isfield(c.myClasf.raw,'fullSgShuf')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShufBal');
    c.myClasf.raw = rmfield(c.myClasf.Zsc,'fullSgShufBal');
end





