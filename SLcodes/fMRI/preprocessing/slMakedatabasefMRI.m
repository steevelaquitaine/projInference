

%slMakedatabasefMRI.m
%
%
% author: steeve laquitaine
%purpose: get and organize data from fMRI experiment
%		  This works for one scan (concatenation of scans or single scan)
%		  Easily have a look at all data (fMRI, behavior and task) collected for a scan
%
%  usage:
% 			o.myGroup 		  = 'MotionComp' 
%        	o.taskNum         = 2;
%        	o.myROIname       = {'rMT'}; 						   %optional
%         	o.myROIpath       = '~/data/mlrAnatDB/s0025/mlrROIs/'   %optional
%           o.myROIpath       = '~/Dropbox/myDropbox/Project_withJustin/data/data/mlrAnatDB/s0025/mlrROIs/'
% 			slMakedatabasefMRI(o)			
%			
%			
%Description:
%
%	- load scan 4D space time BOLD time series (memory and time consuming)
%
%

function [v data o] = slMakedatabasefMRI(o)

%--------
%get data
%--------
fprintf('%s \n','(slMakePredictionsWJM) Gathering data ...')
fprintf('%s \n','(slMakePredictionsWJM) Please set dataPath ...')
o.dataPath = uigetdir(cd,'Pick a fMRI project e.g., ~/data/.../s02520150814...');
cd(o.dataPath)

%The data unit is the smallest., fMRI acquisition vols (TR)
%It might be memory intensive for very small TR but that can be
%alleviated by averaging and storing over periods of interest
%when asked for.

%mrLoadRet
cd(o.dataPath)
mrLoadRet([])
v = getMLRView;

%set scan
v = viewSet(v,'curGroup',o.myGroup);
v = viewSet(v,'curScan',1);

%---------
%fMRI data
%---------
%load 4D Bold data by trials (cells)
data.scan          = loadScan(v);
data.BoldTSeries   = data.scan.data;      %4D space time image (NvoxHeight by NvoxWidth by Nslices by NtimeVols)

%---------------
%Behavioral data
%---------------
%get psychophsyics data from stimfile associated with scan
e = getTaskParameters(data.scan.stimfile{1}.myscreen,data.scan.stimfile{1}.task);
e = e{o.taskNum};  
[~,data.estimatesDeg,~] = SLcart2polar(cell2mat(e.randVars.prodcoor')); %estimates

%---------
%Task data
%---------
data.sStrg          = e.randVars.myRandomCoh;  %coherence
data.s              = e.randVars.myRandomDir;  %direction
data.priorModes     = e.randVars.myModes;      %prior modes
data.pstd           = e.randVars.myStrength;   %prior strength



