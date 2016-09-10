
%usage
%
%     o.sessPath = '/Volumes/Transcend/data/sltaskdotdirfmri05/s02520150814';
%     o.taskNum = 2;
%     o.phaseNum = 1;
%     o.segmentNum = 2;
%     variables = {'myRandomDir','myRandomCoh'}
%     test00(o)


function d = test00(o,variables)

%open session
cd(o.sessPath)
v = mrLoadRet([]);

%get all stimvols and variables
fprintf('\n %s \n \n','(slGetSessStimvol) Displaying stimvols by classes and scans----')
[d,o.taskCond] = getStimvol(v,'_all_','taskNum',o.taskNum,...
    'phaseNum',o.phaseNum,'segmentNum',o.segmentNum);


dbstack
keyboard