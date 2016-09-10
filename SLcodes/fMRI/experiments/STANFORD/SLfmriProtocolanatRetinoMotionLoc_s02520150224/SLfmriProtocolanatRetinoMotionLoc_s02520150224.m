
% SLfmriProtocolanatRetinoMotionLoc_s02520150224.m
% 
%                  $Id: SLfmriProtocolanatRetinoMotionLoc_s02520150209.m 380 2015-01-12 07:00:55Z steeve $
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
%
% 1) CNI 3D anatomy (1 scan, 5 min)
%..................................
% 
% Standard anatomy (GE standard) 0.9 x 0.9 x 0.9 mm voxel resolution, (5 min scan, this session)
% Standard anatomy with acceleration (2x or 3x, 1 - 2 min scan, will be used in subsequent sessions)
%
%
%
% 2) Remove distortions (2 scans, 3 min)
%........................................
% 
%   - Field map (1 - 2 min scan)
%   - PE polar (10 sec)
% 
%
%
% 3) EPI (8 scans, 32 min)
%..........................
%
%8 retinotopy scans of 4 min
%
% expanding ring (10 cycles, 24 s per) - 2X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res
% expanding ring (10 cycles, 24 s per) - 2X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
% expanding ring (10 cycles, 24 s per) - 3X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
% expanding ring (10 cycles, 24 s per) - 3X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
%
%
% contracting ring (10 cycles, 24 s per) - 2X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
% contracting ring (10 cycles, 24 s per) - 2X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
% contracting ring (10 cycles, 24 s per) - 3X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
% contracting ring (10 cycles, 24 s per) - 3X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res.
%
%
%summary
%   -  8 scans of ring stimulus: one for each acceleration sequence X2
%      10 cycles of the ring stimulus for each scan (24 sec period *10 = 4 min)
%      ref: http://gru.stanford.edu/doku.php/mrTools/tutorialsRetinotopy
%   -  2.4 mm x 2.4 mm x 2.4 mm voxel resolution
%   -  4 accelerations sequences: 2 slice accelerations (2x and 3x) X 2 in plane accelerations (1x and 2x)
%
%
%Retinotopy
%----------
%
% Traveling-wave retinotopic mapping experiment and fixation task
% A scan is typically a stimulus ring that revolves over the full 
% circle (a cycle) over 24 sec 10 times (10 cycles lasts 4 min) 
% 
% We run 8 scans that progress in the order:
% - expanding rings (4 min)
% - contracting rings (4 min)
%
% 
% default parameters
% - 10 cycles per scan (e.g., wedges revolve over the full circle 10 times)
% - a cycle's period is 24 sec (scan is thus 4 min)
% - initial half cycle: 1 (add half cycle for response stabilization; will be thrown out)
% - 0.25 duty cycle
% - easyFixTask 1 (size and timing of fixation cross for fixation task)
% - eyeCalb -1 (eye calibration at the start of the scan)
% - dispText: default
% 
% other parameters
% - 16 volumes per cycle (= cycle period./TR  [period = 24 sec; TR = 1.5 sec])
% - 168 volumes per scan (10 cycles * 16 volumes + 8 vols thrown out)
% 
% The code requires a screen parameters file .mglScreenParams.mat and saves 
% a .mat file for each scan in a folder ~data/mglRetinotopy
% The .mat file contains:
% - stimulus parameters
% - screen parameters
% - timing of the volumes acquired by the fMRI scanner
% 
%
%3 - Localizers (current code)
%----------------------------------------
%
%A scan:
%
% - 10 cycles of 16 volumes + 1 cycle to stabilize fMRI
% - 1 cycle alternates 12s expand/contract with 12 sec random coh
% - 176 volumes per scan (88 vols when acceleration x2)
% - TR = 1.57 sec 
% - a scan lasts 4 min 36s                
%
%--------                       
%About TR
%--------
%note: TR may change slightly between 1.5 & >1.5. If the number of volumes 
%      is chosen relative to trial duration to avoid too short TR that
%      do not acquire all the trials.
%      e.g., If I want to collect 8 vols/trials and 1.5 sec <= TR <= 1.71 sec
%      a trial must last 12 sec because 12 sec./8 vols <= TR <= 12./7 vols
%
%--------------
%Testing a scan
%--------------
%note: mglSimulateRun is used only to test the code. It sends a pulse 
%      that triggers volume acquisition. It sometimes causes matlab to crash
%
%      e.g., If you want to simulate how the scanner pulses triggers the 
%      stimulus ('a' means that the code waits for 'a' to start):
%
%            mglSimulateRun(1.57,176,3,'a'); 
%            dotslocSL('randMotion')
%
%
%
%
%COMMAMDS
%--------
%
%Retinotopy: rings
%------------------
%
%Expanding ring (10 cycles, 24 s per) - 2X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR (4 min)
%mglSimulateRun(1.5,168,1);

%mglSimulateRun(1.5,60,1);
% eccArray = [5.5, 7.5, 17, 18, 19, 11, 8, 6,14 , 15, 11, 8, 7.5];
% mglOpen
% stencil = mglProjStencil(eccArray/28.4);
% mglClearScreen(0);
% mglStencilSelect(1);
% mglFillRect(0,0,[25 25],[.5 .5 .5]);
% mglStencilSelect(0);
% mglFlush
%Expanding ring (10 cycles, 24 s per) - 2X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
% mglClose
%ESC

%Expanding ring (10 cycles, 24 s per) - 3X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Expanding ring (10 cycles, 24 s per) - 3X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC


%Contracting ring (10 cycles, 24 s per) - 2X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
%mglSimulateRun(1.5,168,6,'a');
mglRetinotopy('displayName=fMRIproj','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Contracting ring (10 cycles, 24 s per) - 2X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Contracting ring (10 cycles, 24 s per) - 3X slice acc. - 1X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Contracting ring (10 cycles, 24 s per) - 3X slice acc. - 2X in plane accelerations - 2.4 mm x 2.4 mm x 2.4 mm vox.res, 1.5s TR  (4 min)
mglRetinotopy('displayName=fMRIproj','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1'); 
%ESC



%If there's time, Localizer
%-------------------------
%localizer that alternates 100% coh and 0% coh motion dots, 1.57 s TR (4 min 36)
%mglSimulateRun(1.57,176,0,'a'); 
SLfmriLocalizerMotion('randMotion')

%localizer that alternates motion dots amd black screen, 1.57 s TR (4 min 36)
%mglSimulateRun(1.57,176,6,'a'); 
SLfmriLocalizerMotion('black')

