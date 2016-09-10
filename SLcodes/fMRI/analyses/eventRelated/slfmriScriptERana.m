

% slfmriScriptERana.m
%
%author: steeve laquitaine
%  date: 160106
%
%
%%Event-related analyses (stim-locked Bold responses)
%HDR is de-convolved from the scan BOLD time series at each 
%event (e.,g motion onset, a unique HDR for all task conditions 
%calculated for 50 vols, 25 sec with .5s TR) then fitted to the 
%BOLD time series. R^2 assesses the quality of the fit. 
%The error bars are std over \# of repeats of the stimulus 
%(fit residual variance distributed back to each time point 
%in the de-convolved response according to the inverse of the 
%covariance of the design matrix. Thus std is the same at all points).
%analyse occipit(V1,V3A,MT) and parietal(IPS1,2,3)
%task: [fixation(1s)] - [Motion(.3)] - [ISI(5s)] - [Resp(3s)] - [ITI(11.7s)]

%init
o.myVar = 'myRandomCoh';
o.myROIname = {'V1','V3A','MT','IPS'};
o.myGroup = 'Concatenation';
o.myGroupScanNum = 1;
o.scanList = 1:18;
o.myr2Thresh = 0.2;
o.taskNum = 2;
o.phaseNum = 1;
o.segmentNum = 2; %motion
o.hdrlen = 50;    %(in frames, motion + 5 sec)
o.myBase = 's0025_flatR_WM_occipital_Rad90.hdr';
o.myPath = '~/dataold/datafMRI/sltaskdotdirfmri05/s02520150814';
o.myBasePath = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrBaseAnatomies/';
o.myCanonical = '~/dataold/datafMRI/mlrAnatDB/s0025/surfaces/s0025_mprage_pp.nii';
o.myROIpath = '~/dataold/datafMRI/mlrAnatDB/s0025/mlrROIs/';

%run
o = slfMRIerAnaltmp(o);
