
% runAllTaskDirfMRI05.m
%
%  Author : Steeve Laquitaine
%    date : 150813 updated (increased #vols, one loc at start and end, 150924)
% Purpose : Identify areas that represent switching percept (evidence or prior)
% Subject : s025
%   Usage : paste each line. 
% Details : this code requires 20 files
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIOR 225 deg (DAY 1) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ----------- TRAINING (psychophysics outside scanner) ------------
% prior 225 deg, coh: 8 and 16%
% prior: mean: 225 deg - std:80 deg
% 259 trials each
% runs of ~8:30 min
%
sid = 323;
mglSetSID(sid);

%% First, quick practice with 100% coherence (change at ~line 61)
%then real task:
%run 3 times for training to the prior
sltaskDotDirfMRI05_train('VPixx','TRAINparamsPriors80m225Coh12and06',1,3) %12 and 24%
sltaskDotDirfMRI05_train('VPixx','TRAINparamsPriors80m225Coh12and06',1,3) %12 and 24%

% note: s0327: coh was too hard, subject report random so I increased to 12% and 24%
%note: s024: good
%% check data with: cd to stimfile path; load; run "analyses.m"
%with: cd to stimfile path; load; run "analyses.m"
cd(['~/data/sltaskDotDirfMRI05_train/s' num2str(sid)])

%%
analyses({'sub01'},{'StimStrength','Pstd','FeatureSample'},'inpath',...
    'experiment','vonMisesPrior','dataDis','signDistance')





%% ---------------------- fMRI prior 225deg ------------------------
%    #vols : 796 (1 block is 780 ( = 15trials*52vols per trial) + 16 init vols 
%           - 1 vol to make sure the scan stops before the task.
% duration :......
% 15 trials each
%mglSimulateRun(.5,796,5); 
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block01',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block02',1,3)
%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block03',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block04',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block05',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block06',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block07',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block08',1,3)
%% One shorter scan
%    #vols : 380 (364 ( = 7 trials * 52 vols per trial) + 16 init vols - 1 
%            to make sure the scan stops before the task)
% duration : ....
% 7 trials
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m225Coh12and06_block09',1,3)
%% Motion Localizer 
%   #vols : 548
%duration : 4:34 min
% Mtloc('0%',.5) 
            




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIOR 135 deg (DAY 2) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ----------- TRAINING (psychophysics outside scanner) ------------
% prior 135 deg
% mean: 135 deg - std:80 deg
% 259 trials each
% runs of ~8:30 min
%
%First, quick practice with 100% coherence (change at ~line 61)
%then real task:
%run 3 times for training to the prior
sltaskDotDirfMRI05_train('VPixx','TRAINparamsPriors80m135Coh12and06',1,3)%12 and 24%
sltaskDotDirfMRI05_train('VPixx','TRAINparamsPriors80m135Coh12and06',1,3)%12 and 24%

%check data with: cd to stimfile path; load; run "analyses.m"


%% ---------------------- fMRI prior 135deg ------------------------

%% Motion localizer
% dotslocSL('randMotion')
%%
%
%    #vols : 796 (1 block is 780 ( = 15trials*52vols per trial) + 16 init vols 
%           - 1 vol to make sure the scan stops before the task.
% duration :......
% 15 trials each
% mglSimulateRun(.5,796,5) 
%
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block01',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block02',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block03',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block04',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block05',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block06',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block07',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block08',1,3)
%% One shorter scan
%    #vols : 380 (364 ( = 7 trials * 52 vols per trial) + 16 init vols - 1 
%            to make sure the scan stops before the task)
% duration : ....
% 7 trials
% SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsPriors80m135Coh12and06_block09',1,3)
%% Motion Localizer 
%   #vols : 548
%duration : 4:42 min
%Mtloc('0%',.5)            
%569 vols = 553 + 16 initial vols
%mglSimulateRun(.5,569)
%dotslocSL('randMotion')


%% ---------------------- fMRI no prior coh 100 ------------------------

%% Motion localizer
%569 vols = 553 + 16 initial vols
%mglSimulateRun(.5,569)
%dotslocSL('randMotion')

%% Main scans
%    #vols : 796 (1 block is 780 ( = 15trials*52vols per trial) + 16 init vols 
%           - 1 vol to make sure the scan stops before the task.
% duration :......
% 15 trials each
%mglSimulateRun(.5,796,5) 
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block01',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block02',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block03',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block04',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block05',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block06',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block07',1,3)
%%
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block08',1,3)
%% One shorter scan
%    #vols : 276 (260 ( = 5 trials * 52 vols per trial) + 16 init vols
%            to make sure the scan stops before the task)
% duration : 2:18 min
% 5 trials
%SLtaskDotDirfMRI05('fMRIprojFlex','fMRIparamsNoPriorCoh100_block09',1,3)

%% Motion Localizer 
%569 vols = 553 + 16 initial vols
%mglSimulateRun(.5,569)
% dotslocSL('randMotion')


