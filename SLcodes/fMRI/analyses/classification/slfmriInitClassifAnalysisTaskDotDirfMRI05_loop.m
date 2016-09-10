
%slfmriInitClassifAnalysisTaskDotDirfMRI05_loop.m

%% Initialize parameters

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set task and scan parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%o.priormean   = NaN; %for uniform prior
o.priormean     = 225;
%o.priormean   = 135;

o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
o.stckPath = o.rootpath;
o.sessPath{1} = [o.rootpath 's02520150814'];
o.sessPath{2} = [o.rootpath 's02520150923'];
o.sessPath{3} = [o.rootpath 's02520150925'];
o.myGroup = 'Concatenation';
o.myScan = 1;
o.myROInameAll = {'outsideBrain03','V1','V2','V3','V3A','MT','hV4','IPS'};%;{'outsideBrain03','V1'};%,'V2','V3','V3A','hV4','MT','IPS'};
o.taskNum = 2;
o.phaseNum = 1;
o.segmentNum = 2;   
o.anatFileName = 's0025_flatL_WM_occipital_Rad90';
o.anatFilePath = '~/data/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
o.myROIpath = '~/data/datafMRI/mlrAnatDB/s0025/mlrROIs/';

%saved results path
o.classifResultPath  = '/Volumes/DroboBKUP/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif/';

% calculate/add "switch" variable to sessions (once)
[o,behbySess] = slmakeSwitchingVar(o);  



