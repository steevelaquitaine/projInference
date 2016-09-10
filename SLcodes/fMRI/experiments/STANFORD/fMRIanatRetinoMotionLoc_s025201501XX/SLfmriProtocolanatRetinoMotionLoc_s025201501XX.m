
% SLfmriProtocolanatRetinoMotionLoc_s025201501XX.m
% 
%                  $Id: SLfmriProtocolanatRetinoMotionLoc_s025201501XX.m 380 2015-01-12 07:00:55Z steeve $
%                   by: Steeve Laquitaine
%                 date: 15/01/12
%              Purpose: An fMRI experiment with 3D anatomy, retinotopy and motion localizer
% 
%               Usage: run each line for each fMRI scan
%
%Overview
%
% 1 - 3D anatomy (CNI available protocol)
% 2 - Retinotop: Wedges and rings
%     ref: http://gru.stanford.edu/doku.php/mrTools/tutorialsRetinotopy
%   
%           mglDoRetinotopy
%           mglRetinotopy('bars=1')
%
% 3 - Retinotopy: bars
%
% 4 - localizer that alternates 100% coh and 0% coh motion dots
%  
%           SLfmriLocalizerMotion('randMotion')
%
% 5 - localizer that alternates motion dots amd black screen
%  
%           SLfmriLocalizerMotion('black')
%
%note: Requires mgl and mrTools distributions to work.
%
%Details
%
%1 - 3D anatomy (CNI available protocol)
%
%
%2 - Retinotopy (current code)
%------------------------------
% Traveling-wave retinotopic mapping experiment and fixation task
% A scan is typically a stimulus wedges or rings that revolve over the full 
% circle (a cycle) over 24 sec 10 times (10 cycles lasts 4 min) 
% 
% We run 10 scans that progress in the order:
% - CCW wedges (4 min) 
% - CW wedges (4 min) 
% - expanding rings (4 min)
% - contracting rings (4 min)
% - CCW wedges (4 min) 
% - CW wedges (4 min) 
% - bars (4 min) 
% - bars (4 min) 
% - bars (4 min) 
% - bars (4 min) 
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
% - 16 scans
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

%COMMAMDS
%--------
% 1 - 3D anatomy (CNI available protocol)

% 2 - Retinotopy: Wedges and rings
%CCW
%mglSimulateRun(1.5,168,1,'a'); 
% mglSimulateRun(1.5,160,0,'a'); 
mglRetinotopy('wedges=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=1','initialHalfCycle=1');
mglSimulateRun(1.5,160,1,'a'); 
mglRetinotopy('wedges=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%CW
%mglSimulateRun(1.5,168,1,'a'); 
mglRetinotopy('wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Expanding
%mglSimulateRun(1.5,168,1,'a');
mglRetinotopy('rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%Contracting
% mglSimulateRun(1.5,168,6,'a');
mglRetinotopy('rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%CCW
% mglSimulateRun(1.5,168,1,'a');
mglRetinotopy('wedges=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC

%CW
% mglSimulateRun(1.5,168,1,'a'); 
mglRetinotopy('wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%ESC


% 4 - Retinotopy: bars (x4)
%mglSimulateRun(1.5,168,1,'a'); 
mglRetinotopy('bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1')
mglRetinotopy('bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1')
mglRetinotopy('baars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1')
mglRetinotopy('bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1')

% 5 - localizer that alternates 100% coh and 0% coh motion dots
%mglSimulateRun(1.57,176,0,'a'); 
SLfmriLocalizerMotion('randMotion')

% 6- localizer that alternates motion dots amd black screen
%mglSimulateRun(1.57,176,6,'a'); 
SLfmriLocalizerMotion('black')

