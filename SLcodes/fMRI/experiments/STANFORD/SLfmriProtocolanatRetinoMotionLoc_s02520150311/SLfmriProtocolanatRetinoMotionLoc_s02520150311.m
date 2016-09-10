
% SLfmriProtocolanatRetinoMotionLoc_s02520150311.m
% 
%                  $Id: SLfmriProtocolanatRetinoMotionLoc_s02520150311.m 380 2015-01-12 07:00:55Z steeve $
%                   by: Steeve Laquitaine
%                 date: 15/02/24
%              Purpose: An fMRI experiment with 3D anatomy, retinotopy and motion localizer
%               Usage: 
%
%                       run each line for each fMRI scan
%
%note: Requires mgl and mrTools distributions to work.
%
%--------
%Overview
%--------
% Total time: 40 min
%
% 1) CNI 3D  canonical anatomy .9 mm (1 scan, 4-5 min)
%..................................................
% Canonical anatomy (GE standard) 0.9 mm iso voxel resolution
%
%
% 2) test 5 CNI 3D session anatomy 1.2 mm (1 scan, 2 min)
%........................................................
% - sagittal
% - axial
% - coronal
% - parallel to calcarine
% - perpendicular to calcarine
% 
% 3) fMRI retinotopy (8 scans, 32 min)
%......................................
%
% Traveling-wave retinotopic mapping experiment and fixation task
% A scan is a stimulus ring that revolves over the full 
% circle (a cycle) over 24 sec 10 times (10 cycles lasts 4 min) 
% 
% We run 8 expanding ring scans
% 
% default parameters
% - 10 cycles per scan    (e.g., wedges revolve over the full circle 10 times)
% - 24 sec/cycle
% - initial half cycle: 1 (start with an half cycle for response stabilization; will be thrown out)
% - 0.25 duty cycle
% - easyFixTask 1         (size and timing of fixation cross for fixation task)
% - eyeCalb -1            (eye calibration at the start of the scan)
% - dispText: default
% 
% e.g. parameters for TR = 1.5s
% - 16 volumes per cycle (= cycle period./TR  [period = 24 sec; TR = 1.5 sec])
% - 168 volumes per scan (10 cycles * 16 volumes + 8 vols thrown out)
% 
% The code requires .mglScreenParams.mat and saves a .mat file for each scan 
% in a folder ~data/mglRetinotopy_proj with
% - stimulus parameters
% - screen parameters
% - timing of the volumes acquired by the fMRI scanner
%       
%
%
%to test simulate a scan run:
%
%            mglSimulateRun(1.57,176,3,'a'); 
%            mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');
%
%
%
%COMMAMDS
%--------
%
% 1) 1 canonical anatomy : T1w .9mm (CNI software)
% 
% 2) 5 session anatomies:  T1w 1.2 mm (CNI software)
% - sagittal
% - axial
% - coronal
% - parallel to calcarine
% - perpendicular to calcarine
%
%3) fMRI (retinotopy, CNI software and command line)
%
%Expanding ring (4 min each)
% mglSimulateRun(1.4,180,1);
%10 cycles, 24 s/ cycle, MUX2 slice acc. - ARC1 in plane accelerations - 2.4 mm iso vox.res, 1.4s TR, 180 vols acq
%mglSimulateRun(1.4,180,1,'a');

mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');

%10 cycles, 24 s per cy, MUX2 slice acc. - ARC2 in plane accelerations - 2.4 mm iso vox.res, 1.2s TR, 200 vols acq
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');

%10 cycles, 24 s per cy, MUX3 slice acc. - ARC1 in plane accelerations - 2.4 mm vox.res, 1.15s TR, 210 vols acq
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');

%10 cycles, 24 s per, MUX3 slice acc. - ARC2 in plane accelerations - 2.4 mm vox.res, .95s TR, 260 vols acq
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');
mglRetinotopy_proj('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1','projector=1');



