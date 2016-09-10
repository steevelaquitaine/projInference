

function slGetSessionInstances


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