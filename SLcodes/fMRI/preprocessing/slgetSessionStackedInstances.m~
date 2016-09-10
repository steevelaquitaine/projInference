
%slgetSessionStackedInstances.m
%
%author: steeve laquitaine
%  date: 160122
%purpose: calculate instance for each fMRI session and stacked them
%         together by aligning voxels.
%
%usage 
%       [o,si,s,c] = slgetSessionStackedInstances(o)
%
%
%inputs : 
%                o.sessPath: session paths
%                       e.g., 
%                           o.sessPath{1} = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814';
%                           o.sessPath{2} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150923';
%                           o.sessPath{3} = '~/data/datafMRI/sltaskdotdirfmri05/s002520150925';
%              o.myRoiname : e.g., 'MT'
%                  o.nSess : e.g.,3
%                  o.nRois : e.g.,1
%                o.myGroup : e.g.,'Concatenation'
%                 o.myScan : e.g.,1
%                o.myClass : e.g.,'myRandomCoh'
%               o.taskCond : e.g.,{'myRandomCoh=0.06','myRandomCoh=0.12'}
%                o.taskNum : e.g.,2?
%             o.segmentNum : e.g.,2
%o.SavedroiTSeriesfilename : saved roi tseries file path. Much 
%                            faster than re-extracting roi Tseries
%                            with loadROITSeries function. This
%                            load ROI Tseries data that have already 
%                            been loaded with loadROITSeries and 
%                            then saved in a tSeries....mat file
%                           e.g., ~/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/s02520150814/roiTseries/IPS/tSeriesConcatenations02520150814IPS19-Jan-2016.mat
%
%outputs:
%
%   - save roi tSeries in path : [o.stckPath 'slStckAnalyses/']

%calculate instances
function [o,si,s,c] = slgetSessionStackedInstances(o)

fprintf('\n %s \n \n','---- Calculate and stack instances over sessions -------' )

%# of sessions
o.nSess = length(o.sessPath);

%# of Rois
o.nRois = length(o.myROIname);
s = [];

%loop over sessions and get stimvols and roi instances
for si = 1 : o.nSess
    
    %[o,si,s,c] = slGetSessionInstances(o)
    
    %init
    d = [];    
    
    %load this session
    cd(o.sessPath{si})
    mrQuit
    h = mrLoadRet([]);
    
    %session info
    o.si = si;
    o = sldispAnalysisInfo(o);
    
    %set group & scan
    v{si} = h;
    v{si} = viewSet(v{si},'curGroup',o.myGroup);
    v{si} = viewSet(v{si},'curScan',o.myScan);
              
    %---------------- get  stimvols by class -------------
    %checked : all good
    fprintf('\n %s \n \n','(slgetSessionStackedInstances) Stimvols by class and scan ----')
    [d,o.taskCond] = getStimvol(v{si},o.myClass,'taskNum',o.taskNum,...
        'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);         
    
    %check if stimvol field exists
    if isfield(d,'stimvol'); o.stimvolAll = d.stimvol; else o.stimvolAll = d; end
    
    %session basic info
    maxvols = max([o.stimvolAll{:}]);
    fprintf('%s %i \n','(slgetSessionStackedInstances) Max # of vols: ',maxvols)
    fprintf('%s \n','(slgetSessionStackedInstances) Classes are: ')
    numStimVols = nan(length(o.taskCond),1);
    for iClass = 1 : length(o.taskCond)
        numStimVols(iClass) = size(o.stimvolAll{iClass},2);        
        disp(sprintf('Class %i is : %s - %i stimvols',iClass,o.taskCond{iClass},numStimVols(iClass)));        
    end
    
    %When a variable is regressed-out get its stimvols
    if any(strcmp(o.myAnalysis,'regOutVar2'))
        fprintf('%s %s %s %s \n','(slfmriClassify) Variable',o.regOutVar2,'will be regress-out before classification of variable',o.myClass)
        fprintf('%s %s %s %s \n','(slfmriClassify) Getting its stimvols')
        [o.stimvol2,o.taskCond2] = getStimvol(v{si},o.regOutVar2,'taskNum',o.taskNum,'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);
    end        
    
    %------------------- get instances by roi -------------
    [o,s] = slgetROIinstances(o,v,si,s);
    
    %keep original stimvols
    o.stimvolAlls{si} = o.stimvolAll;
    
    %clear stimvolAll field
    o = rmfield(o,'stimvolAll');     
end



%stack over sessions by roi
nVar = length(o.taskCond);
for roii =  1 :  o.nRois
    s = slStackSessInst(o,s,nVar,roii);    
end

%save in directory
[o,s] = slSavedInStrDirectory(o,s);



%display session info
function o = sldispAnalysisInfo(o)

fprintf('\n ----- Session info -------- \n')
fprintf('\n %s %i %s %i \n','(slgetSessionStackedInstances) Session:',o.si,'of',o.nSess)
fprintf('%s %s \n','(slgetSessionStackedInstances) Path:'         ,o.sessPath{o.si})
for i = 1 : o.nRois
    fprintf('%s %s \n','(slgetSessionStackedInstances) ROI: '     ,o.myROIname{i})
end
fprintf('%s %i \n','(slgetSessionStackedInstances) Task: '        ,o.taskNum)
fprintf('%s %i \n','(slgetSessionStackedInstances) Task phase: '  ,o.phaseNum)
fprintf('%s %i \n','(slgetSessionStackedInstances) Task segment: ',o.segmentNum)
o.savedClass = o.myClass;
while iscell(o.savedClass)
    o.savedClass = [o.savedClass{:};];
end
fprintf('%s %s \n','(slgetSessionStackedInstances) Class:'       ,o.savedClass)



%get session instances by roi
function [o,s] = slgetROIinstances(o,v,si,s)

%each roi
for ROIi = 1 : o.nRois
    
    %clean up workspace
    keep = {'s','si','v','o','ROIi','keep','a','i'}; a=whos;
    for i = 1 : length(a)
        if ~sum(strcmp(a(i).name,keep))>0; clear(a(i).name); end
    end
    clear keep; clear a; clear i
    
    %-----------------get roi Tseries -----------------------
    %re-loadROITseries (very long, not recommended except if there is no 
    %already existing roi tseries)
    if ~slIsInput(o.myvarg,'loadSavedROI')
        myROIi = slsaveROItseries(v,o,ROIi,si);        
        %or load saved roi tseries (recommended, much faster)
    elseif slIsInput(o.myvarg,'loadSavedROI')
        [~,myROIi] = slimportROItseries(o,si,ROIi);
    end
    
    %roi info
    fprintf('\n %s %s \n','---- ', o.myROIname{ROIi})
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
    
    %instances at specified time
    if any(strcmp(o.myAnalysis,'accAtTime'))        
        thiss = [];       
        thiss = slgetInstancesAtTime(v{si},o.volsToClass(ROIi,:),...
            o.stimvolAll,o.taskCond,c.myROIinfo,o.sessPath{si});        
        %store session info
        s.session{ROIi}{si}  = thiss.session;
        s.stimvols{ROIi}{si} = thiss.stimvols;
        s.i{ROIi}{si}        = thiss.i;
        s.izsc{ROIi}{si}     = thiss.izsc;               
        s.variable           = thiss.variable;    
    end
end
%update stimvols everywhere (as some stimvols are
%sometimes dropped when getting instances)
o.stimvol{ROIi}{si} = s.stimvols{ROIi}{si};
clear c



%save a session ROI tseries in structured directory
function myROIi = slsaveROItseries(v,o,ROIi,si)

%Must be called for each session.
%
%inputs: 
%   si : session index (e.g., 1s or 2nd or... session)
% ROIi : roi index (e.g., 1s or 2nd or... roi)
%
%Each session has its own roi tseries
%roi with path
% fullroi = [o.myROIpath o.myROIname{ROIi} '.mat'];
%load ROITseries
%This is insanely slow, lots of memory !! make sure edit - preferences - maxBlocksize is low enough
%GET rid of it asap. Prob is er analysis doesn't load and store the
%entire volume but load slice by slice. Thus cannot be used here to get
%the ROI TSeries. The fastest way so far is to save the ROI TSeries once
%run once and reload them instead of loadROITSeries them.
%e.myROI.scanCoords = getROICoordinates(v,myROIname);
myROIi = loadROITSeriesFromPath(v{si},o.myROIname{ROIi},[],[],['roidir=' o.myROIpath]); %load roi tseries (very slow)

%check roi
if isempty(myROIi)
    fprintf('\n')
    slPrintfStr('slsaveROItseries',' Could not find ROI. Path incorrect ?.....................')
    fprintf('\n')
    mrQuit; keyboard
end

%session path and name
[sPath,fn] = fileparts(o.sessPath{si});

%check dir exists for this session
savedDir = [o.stckPath 'slStckAnalyses/' o.myGroup '/' fn '/roiTseries/' o.myROIname{ROIi}];
if ~isdir(savedDir)
    mkdir(savedDir)   
end

%save ROI tseries data
cd(savedDir)
save(['tSeries' o.myGroup fn o.myROIname{ROIi} date],'myROIi')
cd ..

%get instances at specified time
function s = slgetInstancesAtTime(v,volsToClass,stimvol,variable,myROIinfo,sessPath)

%block length
c.blockLen = volsToClass(2) - volsToClass(1) + 1;

%instance vols info
fprintf('\n %s \n','------- Averaging instances over vols -------')
fprintf('%s %i \n','# of vols:',c.blockLen)
fprintf('%s %i %s %i %s \n','from event +',volsToClass(1),'to',volsToClass(2), 'vols')

%lag to event
%first vols
c.startLag = volsToClass(1);
c.startLag = ['startLag=' num2str(c.startLag)];
%# of avergaed vols
c.blockLen = ['blockLen=' num2str(c.blockLen)];

%average vols after stim. onset
%for all voxels (n=inf)
Inst = getInstances(v,myROIinfo,stimvol,'n=inf',c.startLag,c.blockLen);
Inst = Inst{1};

%raw instances
s.i = Inst.classify.instances;  

%z-scored instances (over trials)
s.izsc = preprocessInstances(s.i,'zscore=1');

%variable by which instances are sorted
s.variable = variable;

%session name
[~,s.session] = fileparts(sessPath);

%instance stimvols
s.stimvols = Inst.classify.instanceVol;
s.descriptions = {'for each roi: session: session name; simvols: sessions stimvols; i: instances; izsc: zscored-instances'};


%stack instances across sessions
function s = slStackSessInst(o,s,nVar,roii)

%loop over variables and stack instances over session
for vi = 1 : nVar    
    
    %init
    s.iStacked{roii}{vi}  = [];
    s.izStacked{roii}{vi} = [];
    s.stimvolsStacked{roii}{vi} = [];
    
    %stack instances over sessions
    for si = 1 :  length(o.sessPath)
        
        %init
        is = [];
        isz = [];
        
        %warning missing variables
        if length(s.i{roii}{si}) < nVar
            fprintf('%s \n','!!!! There are missing variables !!! Fix that......')
            keyboard
        end        
        %deal with missing classes
%         if vi > length(s.i{roii}{si})
%             s.i{roii}{si}{vi}={};
%         end
%         if vi > length(s.izsc{roii}{si})
%             s.izsc{roii}{si}{vi}={};
%             isz=[];
%         end
        is = s.i{roii}{si}{vi};
        isz = s.izsc{roii}{si}{vi};
%         if isempty(s.i{roii}{si}{vi})
%             is =[];
%         end
%         if isempty(s.izsc{roii}{si}{vi})
%             isz =[];
%         end
        
        %Instances (col) by vox (rows) stacked across sessions
        %sorted by classes (cells)
        s.iStacked{roii}{vi} = [s.iStacked{roii}{vi}   ; is]; %raw
        s.izStacked{roii}{vi} = [s.izStacked{roii}{vi} ; isz];%z-scored
                
        %stack stimvols too
        s.stimvolsStacked{roii}{vi} = [s.stimvolsStacked{roii}{vi} s.stimvols{roii}{si}{vi}];
    end
end

%save stacked instances in structured directory
function [o,s] = slSavedInStrDirectory(o,s)

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
save(filename,'o','s')
