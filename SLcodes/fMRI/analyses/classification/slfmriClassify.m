
%slfmriClassify.m
%
%    Author: steeve laquitaine
%      date: 151225
%   purpose: decoding various features with various options from ROI 
%            voxels Bold pattern
%
%inputs:
%
%       'CalcInstances'(or 'loadSavedInstances')
%       'slReSampInst': re-sample instances when few only
%                       set random values for missing classes. Do not run
%                       when classes are missing !!
%               'test': just testing: save  in a test folder
%
%Loading:
%         'loadSavedROI' : load existing roi TSeries file loaded with
%                          loadROITSeries which is time expensive
%
%cross-validation:
%          'leaveOneOut' : train on all instance-1 and test on left over
%              'kFoldCV' : divide into 10 folds train on 9 folds test on 1
%                          left over
%
%
%Chance level:
%        'chanceByPerm' : calculate chance label from N-permuted dataset
%            'nPerm=100': to get chance level form 100 dataset of permuted
%                         classes
%     'permutationUnBal': permute instances classes
%       'permutationBal': equalize number of instances and permute classes
%
%denoising:
%       'r2cutoff',0.03: get rid of voxels with low r2 in an event-related
%                        analysis. You have to indicate before the name of 
%                        the er file that should be loaded (e.g., o.erAnalFile =
%                        'erAnalMotionByAll.mat')
%
%
%outputs:

function slfmriClassify(o,myClass,varargin)

%duration
t0 = tic;

%getArgs
o.nPerm = [];
getArgs(varargin,{'nPerm','accuracyAtTime','CalcInstances',...
    'slReSampInst','loadSavedROI','kFoldCV'});

%case permutation
if ~ieNotDefined('nPerm')
    o.nPerm = nPerm;
end

%get classes
o.myClass = myClass;

%get analyses
o.myAnalysis = [];


o = slgetAnalysis(o,varargin);

%stack instances over sessions
[o,si,s,c] = slStackSessionInstances(o,varargin{:});

%classify by roi
nVar = length(o.taskCond);
for roii = 1 : o.nRois
    fprintf('%s \n',['(leaveOneOut)' o.myROIname{roii}] )
    
    %path
    cd([o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]])
    
    %init
    %[AligniThisRoi,AlignizThisRoi,s] = slStackSessInst(o,s,nVar,roii);    
    AligniThisRoi = c.Aligni{roii};
    AlignizThisRoi = c.Aligniz{roii};

    %---------------------------- for debugging ----------------
    %simulate instances for test
    %c=simInstances({'lV1'},nVar,[],repmat(25,nVar,1),'type=a2Dclusters');
    %AligniThisRoi=repmat({c{1}.classify.instances},o.nRois,1);
    %AlignizThisRoi=preprocessInstances(AlignizThisRoi{1},'zscore=1');
    %AlignizThisRoi=repmat({AlignizThisRoi},o.nRois,1);
    %---------------------------- for debugging ----------------
    
    %---------------------------- for debugging ----------------
    %increase small dataset by resampling
    %size of dataset.just for test: might break
    %the training/test isolation requirement.
    %also I set to rand the classes when instances
    %that have no instances
    %this produces high %correct
    if slIsInput(varargin,'slReSampInst')
        [AligniThisRoi,AlignizThisRoi] = slReSampInst(AligniThisRoi,AlignizThisRoi);
    end
    %---------------------------- for debugging ----------------
    
    fprintf('\n')
    fprintf('%s \n \n','------ Classification (actual data) --------')
    c.myClasf.raw.fullSg = [];
    c.myClasf.Zsc.fullSg = [];
    %case leaveoneout
    if slIsInput(varargin,'leaveOneOut')        
        c.myClasf.raw.fullSg = leaveOneOut(AligniThisRoi{roii},['type=' o.type],o.dataBalancing);%raw
        c.myClasf.Zsc.fullSg = leaveOneOut(AlignizThisRoi{roii},['type=' o.type],o.dataBalancing);%zs
    end
    
    %case k-fold
    if slIsInput(varargin,'kFoldCV')
        fprintf('%s \n','kFoldCV')
        fprintf('%s \n \n','-----------')        
        c.myClasf.raw.fullSg = kFold(AligniThisRoi{roii},'numFolds=9',['type=' o.type],o.dataBalancing);%raw
        disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',o.type,c.myClasf.raw.fullSg.correct*100));
        
        c.myClasf.Zsc.fullSg = kFold(AlignizThisRoi{roii},'numFolds=9',['type=' o.type],o.dataBalancing);%zs
        disp(sprintf('(kFold) %s classifier produced %0.2f%% correct',o.type,c.myClasf.Zsc.fullSg.correct*100));
    end 
    fprintf('\n')
    fprintf('%s \n \n','---- (Classification at chance) permuted data ---')
    %classify permuted data N times (typically N=100)
    %36 min (50 perms,1 ROI) (Markus Ojala, see Permutation 
    %test for stuyding classifer performance, Journal of 
    %Machine learning, 2010) 
    o.instances = AligniThisRoi{roii};
    o.Zscinstances = AlignizThisRoi{roii};
    [c,o] = slGetClassifChance(o,c);

    %fprintf('%s \n','------ Classification (balanced data) ---------')
    %uncomment. accuracy doesn't change much when we balance the dataset
    %balance data (removing) and classify
    c.myClasf.raw.fullSgShufBal = [];
    c.myClasf.Zsc.fullSgShufBal = [];
    %c.myClasf.raw.fullSgShufBal = leaveOneOut(AligniThisRoi{roii},'balancByRemovI=1');%raw
    %c.myClasf.Zsc.fullSgShufBal = leaveOneOut(AlignizThisRoi{roii},'balancByRemovI=1');%zscore
    
    %save results classification
    [o,c,s,roii] = slSaveClassiRes(o,c,s,roii);
    
    %duration
    duration=toc(t0);
    fprintf('%s %.2f \n','One more roi --------------',duration)
end

%Go to analysis path
% o.analysisPath = [o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.savedClass '/' [o.myAnalysis{:}]];
% cd(o.analysisPath)


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
%if testing 
if sum(strcmp(varargin,'test'))
    o.test=1;
else
    o.test=0;
end
%--------- cross-validation -----------
%leave-one instance out
if sum(strcmp(varargin,'leaveOneOut'))
    o.myAnalysis{end+1} = 'leaveOneOut';
    o.leaveOneOut = 1;
    o.crossVal = 'leaveOneOut';
else
    o.leaveOneOut = 0;
end
%k-Fold
if sum(strcmp(varargin,'kFoldCV'))
    o.myAnalysis{end+1} = 'kFoldCV';
    o.kFoldCV = 1;
    o.crossVal = 'kFoldCV';
else
    o.kFoldCV = 0;
end
%leave-one k mean instance out
%group instances into K average instances by class and leave-one-out
%classification. The idea is to average out past trial Bold confound.
if sum(strcmp(varargin,'slleaveOneOfkmeanInstanceOut'))
    o.myAnalysis{end+1} = 'slleaveOneOfkmeanInstanceOut';
    pos = find(strcmp(varargin,'slleaveOneOfkmeanInstanceOut')) + 1;
    o.kInstances = varargin{pos};
    o.crossVal = 'slleaveOneOfkmeanInstanceOut';
end
%--------------classifier ------------
%check first
if ~any(strcmp(varargin,{'fisher'})) && ~any(strcmp(varargin,{'svm'}))
    fprintf('%s \n','(slgetAnalysis) Please input classifier: "fisher" or "svm"')
    dbstack
    keyboard
end
if sum(strcmp(varargin,'svm'))
    o.myAnalysis{end+1} = 'svm';
    o.type = 'svm';
end
if sum(strcmp(varargin,'fisher'))
    o.myAnalysis{end+1} = 'fisher';
    o.type = 'fisher';
end
%---------- denoising -------------------
if sum(strcmp(varargin,'r2cutoff'))
    o.myAnalysis{end+1} = 'r2cutoff';
    pos = find(strcmp(varargin,'r2cutoff')) + 1;
    o.r2thres = varargin{pos};    
end
%regress out a variable (not yet correct)
if any(strcmp(varargin,'regOutVar2'))
    o.myAnalysis{end+1} = 'regOutVar2';
    pos = find(strcmp(varargin,'regOutVar2')) + 1;
    o.regOutVar2 = varargin{pos};
end



%----------balancing dataset -------------
if any(strcmp(varargin,'balancByRemovI=1'))
    o.myAnalysis{end+1} = 'balancByRemovI';    
    o.dataBalancing = 'balancByRemovI=1';
elseif any(strcmp(varargin,'balancByRemovI=0'))
    o.dataBalancing = 'balancByRemovI=0';
end
if any(strcmp(varargin,'balancByBootSt=1'))
    o.myAnalysis{end+1} = 'balancByBootSt';    
    o.dataBalancing = 'balancByBootSt=1';
elseif any(strcmp(varargin,'balancByBootSt=0'))
    o.dataBalancing = 'balancByBootSt=0';
end

%---------- chance level ----------------
%calculate chance by permutation
if sum(strcmp(varargin,'chanceByPerm'))
    o.myAnalysis{end+1} = '_chanceByPerm';
    o.chanceByPerm = 1;
else
    o.chanceByPerm = 0;
end
%balancing
if any(strcmp(varargin,'permutationUnBal=1'))
    o.myAnalysis{end+1} = 'permutationUnBal';    
    o.chanceBalancing = 'permutationUnBal=1';
end
if any(strcmp(varargin,'permutationBal=1'))
    o.myAnalysis{end+1} = 'permutationBal';    
    o.chanceBalancing = 'permutationBal=1';
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
fprintf('%s %s \n','(slfmriClassify) Class:',o.savedClass)

%stack instances over sessions
function [o,si,s,c] = slStackSessionInstances(o,varargin)

fprintf('%s \n','---------------------------------------------------')
fprintf('%s \n','---- Calculate instances and stack over sessions ---')
fprintf('%s \n','---------------------------------------------------')
if slIsInput(varargin,'CalcInstances')
    %get stacked instances
    [o,si,s,c] = slCalcStackedInstances(o);
    %backup
    [o,c,s] = slSavedInStrDirectory(o,c,s);
    %or load
elseif slIsInput(varargin,'loadSavedInstances')
    %convert cell class to string
    o.savedClass = o.myClass;
    while iscell(o.savedClass)
        o.savedClass = [o.savedClass{:};];
    end
    cd([o.stckPath 'slStckAnalyses/' o.myGroup 'classif' o.savedClass [o.myAnalysis{:}]])
    %check file exists
    l = dir(['ClassifStckSessData' o.myClass '_allRois_*']);
    if isempty(l)
        slPrintfStr('slStackSessionInstances',...
            [' ClassifStckSessData' o.savedClass '_allRois_..not found...'])
        dbstack
        keyboard
    else
        load(l.name)
    end
else
    slPrintfStr('slStackSessionInstances','Please input either ')
    slPrintfStr('slStackSessionInstances','- loadSavedInstances : to load existing instances')
    slPrintfStr('slStackSessionInstances','- CalcInstances : to calculate instances')
    dbstack
    keyboard
end

%calculate instances
function [o,si,s,c] = slCalcStackedInstances(o)

%inputs:
% o.sessPath
% o.myRoiname
% o.nSess
% o.nRois
% o.myGroup
% o.myScan
% o.taskCond
% o.myClass
% o.taskNum
% o.segmentNum
% o.stimvol

%mrLoadRet
tic

%get instance by session
o.nSess = length(o.sessPath);
o.nRois = length(o.myROIname);
s = [];

for si = 1 : o.nSess
        
    %go to data
    cd(o.sessPath{si})
    
    %have to load each new session
    h = mrLoadRet([]);
    fprintf('\n %s \n \n', ['(slCalcStackedInstances) Session:' o.sessPath{si} ' ------'])
    
    %disp info
    o = sldispAnalysisInfo(o);
    
    %set group & scan
    v{si} = h;
    v{si} = viewSet(v{si},'curGroup',o.myGroup);
    v{si} = viewSet(v{si},'curScan',o.myScan);
    %refreshMLRDisplay(v{si}); %update display
    
    %get this session stimvols
    o = slGetSessStimvol(o,v,si);
    
    %get this session roi instances
    [o,v,si,s] = slgetROIinstances(o,v,si,s);    
    mrQuit
end

%stack over sessions by roi
nVar = length(o.taskCond);
for roii =  1 :  o.nRois
    [c.Aligni{roii},c.Aligniz{roii}] = slStackSessInst(o,s,nVar,roii);    
end

%get session stimvol by class
function o = slGetSessStimvol(o,v,si)

%get stimvols by class
fprintf('\n %s \n \n','(slGetSessStimvol) Displaying stimvols by classes and scans----')
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
    disp(sprintf('(slGetSessStimvol) Class %i is : %s',i,o.taskCond{i}));
    numStimVols(i) = size(o.stimvol{i},2);
    disp(sprintf('(slGetSessStimvol) Class %i has %i stimvols',i,numStimVols(i)));
end


%if variable is regressed out get its stimvol
if any(strcmp(o.myAnalysis,'regOutVar2'))
    fprintf('%s %s %s %s \n','(slfmriClassify) Variable',o.regOutVar2,'will be regress-out before classification of variable',o.myClass)
    fprintf('%s %s %s %s \n','(slfmriClassify) Getting its stimvols')
    [o.stimvol2,o.taskCond2] = getStimvol(v{si},o.regOutVar2,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
end

%get session instances by roi
function [o,v,si,s] = slgetROIinstances(o,v,si,s)

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
        %loading er anal
        cd([o.sessPath{si} '/' o.myGroup]); 
        v{si} = loadAnalysis(v{si},['/erAnal/' o.erAnalFile]);
        %keep voxels with high r2 (removes noisy features)
        c.myROIinfo = slSortROIvoxels(myROIi,r2cutoff);
    end
    
    %instances at a time
    if any(strcmp(o.myAnalysis,'accAtTime'))
        [s,v,o,c] = slgetInstancesAtTime(v,o,c,ROIi,si,s);
    end
    c = rmfield(c,'myROIinfo');    
end
clear c

%save a session ROI tseries in structured directory
function myROIi = slsaveROItseries(v,o,ROIi,si)
%Must be called for each session.
%Each session has its own roi tseries
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
    slPrintfStr('slsaveROItseries',' Could not find ROI. Change path ?')
    mrQuit; keyboard
end
%session
[sPath,fn] = fileparts(o.sessPath{si});
%check dir exists for this session
savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/' fn '/roiTseries/' o.myROIname{ROIi}];
if ~isdir(savedDir)
    mkdir(savedDir)   
end
%backup
cd(savedDir)
save(['tSeries' o.myGroup fn o.myROIname{ROIi} date],'myROIi')
cd ..

%load ROI tseries in workspace
function [thistSeries,myROIi] = slimportROItseries(o,si,ROIi)

%session
[sPath,fn] = fileparts(o.sessPath{si});
%move to roi path
savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/' fn '/roiTseries/' o.myROIname{ROIi}];
cd(savedDir)
%check that there only is one file
f = dir([savedDir '/tSeries*']);
%get most recent file
lastfile = slGetMostRecentFileIndir;
if length([f.isdir])>1;
    %move older files to archive
    mkdir archive    
    for i = 1 : length(f)
        if ~strcmp(f(i).name,lastfile)
            movefile(f(i).name,[savedDir '/archive'])
        end
    end
    fprintf('%s \n','(slImportROItseries) Getting most recent roi tseries file.')
end
%info dir
thistSeries = dir([savedDir '/' lastfile]);
%import roi data
myROIi = importdata([savedDir '/' thistSeries.name]);
cd(o.sessPath{si})

%get instances at time
function [s,v,o,c] = slgetInstancesAtTime(v,o,c,ROIi,si,s)
%get instances
fprintf('\n %s \n','------- Averaging instances at period -------')
%block length
c.blockLen = o.volsToClass(ROIi,2) - o.volsToClass(ROIi,1) + 1;
%lag to event
c.startLag = o.volsToClass(ROIi,1);
c.startLag = ['startLag=' num2str(c.startLag)];
c.blockLen = ['blockLen=' num2str(c.blockLen)];
%average vols after stim. onset
%and save for each session si
c.myROIinfo1     = getInstances(v{si},c.myROIinfo,o.stimvol,'n=inf',c.startLag,c.blockLen);
c.myROIinfo1{si} = c.myROIinfo1{1};
c.i{si}          = c.myROIinfo1{si}.classify.instances;  
[c.izsc{si},c.pStgs] = preprocessInstances(c.i{si},'zscore=1');
%store session stimvols and instances
%session name
[~,s.session{ROIi}{si}]= fileparts(o.sessPath{si});
%stimvols
s.stimvols{ROIi}{si}   = c.myROIinfo1{si}.classify.instanceVol;
%raw instances
s.i{ROIi}{si} = c.i{si};
%zscored instances
s.izsc{ROIi}{si} = c.izsc{si};
s.descriptions = {'for each roi: session: session name; simvols: sessions stimvols; i: instances; izsc: zscored-instances'};
%backup
o.instBySes_ROI{ROIi}{si} = s.i{ROIi}{si};
o.zscinstBySes_ROI{ROIi}{si} = s.izsc{ROIi}{si};

%stack instances across sessions
function [AligniThisRoi,AlignizThisRoi,s] = slStackSessInst(o,s,nVar,roii)

%loop over variables
for vi = 1 : nVar
    %collect instances each session
    AligniThisRoi{roii}{vi} = [];
    AlignizThisRoi{roii}{vi} = [];
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
        AligniThisRoi{roii}{vi} = [AligniThisRoi{roii}{vi} ; is]; %raw
        AlignizThisRoi{roii}{vi} = [AlignizThisRoi{roii}{vi} ; isz];%z-scored
    end
end

%resample instances to increase
%sample size
function [AligniThisRoi,AligniThisRoiz] = slReSampInst(Aligni,AlignizThisRoi)

%get # vox
nClasses = length(AligniThisRoi{:});
nvox=[];
for cl = 1 : nClasses
    nvox = max([nvox size(AligniThisRoi{:}{cl})]);
end

nReSamp = 33;
for cl = 1 : nClasses
    nI = size(AligniThisRoi{:}{cl},1);
    %if no data
    if isempty(AligniThisRoi{:}{cl})
        AligniThisRoi{:}{cl} = rand(20,nvox);%raw
        AlignizThisRoi{:}{cl} = rand(20,nvox);%zscore
    else
        nReSampix = randsample(1:nI,nReSamp,'true');
        AligniThisRoi{:}{cl} = AligniThisRoi{:}{cl}(nReSampix,:);%raw
        AlignizThisRoi{:}{cl} = AlignizThisRoi{:}{cl}(nReSampix,:);%zscore
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
%if not testsing save data in structured dir
if o.test~=1
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
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '_' date '.mat'],'chance')
    clear chance;

elseif o.test==1
    %otherwise save in test folder
    %backup
    if ~isdir(o.myROIname{roii}); mkdir(o.myROIname{roii});end
    cd(o.myROIname{roii})
    mkdir test; cd test    
    o.myAnalPath2{roii} = pwd;
    %archive existing file
    list = dir('*.mat');
    if ~isempty(list)
        if isempty(dir('Archive')); mkdir Archive; end
        for i=1:length(list); movefile(list(i).name,'Archive/'); end
    end
    save(['ClassifStckSess' o.savedClass '_' o.myROIname{roii} '_' date '.mat'],'o','c','s')
    %save chance accuracy results
    %in a different file
    chance.raw = c.myClasf.raw.fullSgShuf;
    chance.Zsc = c.myClasf.Zsc.fullSgShuf;   
    save(['ClassifStckSessChance' o.savedClass '_' o.myROIname{roii} '_' date '.mat'],'chance')
    clear chance;
end

%cleanup
if isfield(c.myClasf.raw,'fullSgShuf')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShuf');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShuf');
end
if isfield(c.myClasf.raw,'fullSgShufBal')
    c.myClasf.raw = rmfield(c.myClasf.raw,'fullSgShufBal');
    c.myClasf.Zsc = rmfield(c.myClasf.Zsc,'fullSgShufBal');
end

%calculate chance level
function [c,o] = slGetClassifChance(o,c)

%inputs
%   o.crossVal (e.g., 'leaveOneOut' or 'kFold')
%   o.nPerm    (e.g., 2)

%init
c.myClasf.raw.fullSgShuf = [];
c.myClasf.Zsc.fullSgShuf = [];

if o.chanceByPerm==1
    %case leaveOneOut
    if strcmp(o.crossVal,'leaveOneOut')
        if ~isempty(o.nPerm)
            fprintf('%s %i %s \n','(slfmriClassify) Running ',o.nPerm,'permutations...')
            nPerm2 = ['nPerm=' num2str(o.nPerm)];
            type = ['type=' o.type];
            balancing = o.chanceBalancing;
            c.myClasf.raw.fullSgShuf = leaveOneOutNpermut(o.instances,nPerm2,type,balancing);
            c.myClasf.Zsc.fullSgShuf = leaveOneOutNpermut(o.Zscinstances,nPerm2,type,balancing);
        end
    end
    %case ~10-Fold
    if strcmp(o.crossVal,'kFoldCV')                
        if ~isempty(o.nPerm)
            fprintf('%s %i %s \n','(slfmriClassify) Running',o.nPerm,'permutations to get chance.')
            nPerm2 = ['nPerm=' num2str(o.nPerm)];
            type = ['type=' o.type];
            balancing = o.chanceBalancing;
            c.myClasf.raw.fullSgShuf = kFoldNpermut(o.instances,nPerm2,'numFolds=9',type,balancing);
            c.myClasf.Zsc.fullSgShuf = kFoldNpermut(o.Zscinstances,nPerm2,'numFolds=9',type,balancing);
        end
    end
end



