
%slfmriInitClassifAnalysisTaskDotDirfMRI05.m
%
%
%
% author: steeve laquitaine
%purpose: initialize the data info for the classification analyses

clear

%Set task and scan parameters

%motion direction prior parameter
%o.priormean   = NaN; %for uniform prior
o.priormean    = 225;
%o.priormean   = 135;

%enter path of the sessions which data you want to stack
o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';
o.stckPath = o.rootpath;
o.sessPath{1} = [o.rootpath 's02520150814'];
o.sessPath{2} = [o.rootpath 's02520150923'];
o.sessPath{3} = [o.rootpath 's02520150925'];

%scan info
o.myGroup     = 'Concatenation';
o.myScan      = 1;

%Rois
o.myROIname   = {'MT'};%{'outsideBrain02'};%{'V1'};%,'V2','V3','V3A','hV4','MT','V7','IPS','parietalBA39','FEF','vlPFC','dlPFC','vmPFC','OFC','aPFC'};

%task info
o.taskNum     = 2;
o.phaseNum    = 1;
o.segmentNum  = 2;   

%display info
o.anatFileName = 's0025_flatL_WM_occipital_Rad90';
o.anatFilePath = '~/data/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
o.myROIpath    = '~/data/datafMRI/mlrAnatDB/s0025/mlrROIs/';
%o.erAnalFile = 'erAnalMotionByAll.mat';

%saved results path
o.classifResultPath = '/Volumes/DroboBKUP/data/datafMRI/sltaskdotdirfmri05/slStckAnalyses/Concatenation/classif';

% calculate/add "switch" variable to sessions (once)
% [o,behbySess] = slmakeSwitchingVar(o);  



