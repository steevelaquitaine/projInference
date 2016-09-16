

%slfmriGetMotLocr2.m
%
%author: steeve laquitaine
%
% usage: 
%
%       sesspath = '~/data/datafMRI/sltaskdotdirfmri05/s02520150814/'		
%       group = 'Averages';
%       roipath = '~/data/datafMRI/mlrAnatDB/s0025/mlrROIs/V1.mat';
%       [r2,o] = slfmriGetMotLocr2(roipath,sesspath,group);

function [r2,o] = slfmriGetMotLocr2(roipath,sesspath,myGroup,varargin)

%move to session
cd(sesspath)

%load session
v = mrLoadRet([]);

%quick load -reload ROI
[o.myROIpath,o.myROIname] = fileparts(roipath);
o.sessPath =  sesspath;
o.myGroup  = myGroup;
o.stckPath = sesspath;
myROI = loadOrReloadROITseries(v,1,1,o,varargin{:});

%converting the scan coordinates of the ROI to linear coordinates
scanDims = viewGet(getMLRView,'scanDims');
myROI.linearScanCoords = sub2ind(scanDims,myROI.scanCoords(1,:),myROI.scanCoords(2,:),myROI.scanCoords(3,:));

%get r2 map
%get this group's scans and set the view to the 
%desired group and scans
nScans = viewGet(v,'nScans',o.myGroup);
v = viewSet(v,'curGroup',o.myGroup,'curScan',nScans);

%load coherence analysis if exists and get the r2
%r2 describes how well voxel response is fit by a sinewave with 
%the same period as the onset/offset of a localizer stimulus (e.g., motion)
%which indicates how correlated the voxel response is to the stimulus
v = loadAnalysis(v,'corAnal/MotionLocalizer.mat'); %load analysis
overlayNames = viewGet(v,'overlayNames');       
OverlayNum = find(strcmp('co',overlayNames));
viewSet(v,'curOverlay',OverlayNum);
r2 = viewGet(getMLRView,'overlayData',1);

%select this ROI r2
r2 = r2(myROI.linearScanCoords);

