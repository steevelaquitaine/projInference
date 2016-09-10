
% runAllTaskDirfMRI05.m
%
%  Author : Steeve Laquitaine
%    date : 150813 updated (increased #vols, one loc at start and end, 150924)
% Purpose : Identify areas that represent switching percept (evidence or prior)
% Subject : s025
%   Usage : paste each line. 
% Details : this code requires 2 files
% 
%           - fMRIParamsPrior80Coh12and06.mat (parameters from SLinitRunExpUniPriorfMRI.m)
% 
%           - slTaskDotDirfMRI05.m (task)
% 
%
%% Motion Localizer
%#vols : 560 (544 + 16 init vols)
%duration : 278 (4:39 min)
Mtloc('0%',.5)                                                                                             
%
%
%% Main Scans
%------------
%    #vols : 650 (1 block, 630 + 16 init vols)
% duration : 325 sec (5:24 min)                 
%
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block01',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block02',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block03',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block04',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block05',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block06',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block07',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block08',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block09',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block10',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block11',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block12',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block13',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block14',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block15',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block16',1)
%%
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block17',1)

%---------------
%% One more scan
%---------------
%    #vols : 190 (174 + 16 init vols)
% duration : 95 sec (1:34 min)
slTaskDotDirfMRI05('fMRIprojFlex','fMRIParamsPrior80Coh12and06_block18',1)


%% Motion Localizer
%#vols : 560
%duration : 278 (4:39 min)
Mtloc('0%',.5)                     

%%
%#vols : 560
%duration : 278 (4:39 min)
Mtloc('0%',.5)                                                                                             
