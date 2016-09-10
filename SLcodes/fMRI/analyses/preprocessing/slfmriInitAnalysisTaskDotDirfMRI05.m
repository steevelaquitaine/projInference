

%slfmriInitAnalysisTaskDotDirfMRI05.m
%
%
%     author: steeve laquitaine
%    purpose: initialize parameters for analysis of sltaskdotdirfmri05 project
%             make sure workspace always clear !
%
%Description: one roi at a time

%initialize with data info
clear; mrQuit
o.rootpath = '~/data/datafMRI/sltaskdotdirfmri05/';

%%%%%%%%%%%%%%%%%%
% choose subject %
%%%%%%%%%%%%%%%%%%
sid = 25;
sid = num2str(sid); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a/many session/s %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
o.sessPath{1}  = [o.rootpath 's02520150814'];  %s0025 - prior225 - coh6_12
o.sessPath{2}  = [o.rootpath 's02520150923']; %s0025 - prior225
o.sessPath{3}  = [o.rootpath 's02520150925']; %s0025 - prior225
%o.sessPath{1}  = [o.rootpath 's02520160614']; %s0025 - prior135
%o.sessPath{1}  = [o.rootpath 's032720160326']; %s0327 - prior225 - coh8_16
%o.sessPath{1}  = [o.rootpath 's032720160328']; %s0327 - prior135 - coh8_16
%o.sessPath{1}  = [o.rootpath 's02420160427']; %s0024 - prior225 - coh8_16
%o.sessPath{1}  = [o.rootpath 's02420160520'];  %s0024  - prior135 - coh8_16
%o.sessPath{1}  = [o.rootpath 's02520160524'];  %s0025 - priorUnif - coh100

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set task and scan parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prior
%o.priormean   = NaN; %for uniform prior
o.priormean    = 225;
%o.priormean   = 135;

%coherence
o.coh2plot   = 0.06;
% o.coh2plot   = 0.08;
% o.coh2plot   = 0.12;
% o.coh2plot   = 0.16;
% o.coh2plot   = 1;

o.motDir        = 85; %(note that: vox select plot does not depend on this input)
o.myGroup      = 'Concatenation';
o.myScan       = 1;
o.taskNum      = 2;
o.phaseNum     = 1;
o.segmentNum   = 2;   
o.myROIpath    = ['~/data/datafMRI/mlrAnatDB/s0' sid '/mlrROIs/'];
% o.anatFileName = ['s00' sid '_flatL_WM_occipital_Rad90'];
%o.anatFileName = 's0327_right_occipital_r90';
%o.anatFileName = 's0024_right_occipital_r90';
%o.anatFilePath = ['~/data/datafMRI/mlrAnatDB/s00' sid '/mlrBaseAnatomies/'];
%o.erAnalFile   = 'erAnalMotionByAll.mat'; %if r2 cutoff of voxels

%initialize analysis
o.myROIname   = {'outsideBrain03'};%{'V1'}%{'audit_Control','V1','V2','V3','V3A','hV4','MT','V7','IPS','parietalBA39','FEF','vlPFC','dlPFC','vmPFC','OFC','aPFC'};
o.filename    = ['ClassifStckSessmyRandomCoh_' o.myROIname{:} '.mat'];
o.myVar       = 'myRandomCoh'; 
varargin      = {'accuracyAtTime',[7 14],'loadSavedROI','CalcInstances','leaveOneOut','fisher','balancByRemovI=1'};
%varargin     = {'accuracyAtTime',[7 14],'CalcInstances','leaveOneOut','fisher','balancByRemovI=1'};
o             = slfMRIgetAnalyses(o,varargin);
o.stckPath    = o.rootpath;
o.dataPath    = [o.stckPath 'slStckAnalyses/' o.myGroup '/classif/' o.myVar '/' [o.myAnalysis{:}] '/' o.myROIname{:} '/'];


