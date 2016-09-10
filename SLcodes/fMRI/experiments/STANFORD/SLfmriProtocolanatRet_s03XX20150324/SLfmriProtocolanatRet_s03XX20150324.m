
% SLfmriProtocolanatRetinoMotionLoc_s03xx20150324.m
% 
%                  $Id: SLfmriProtocolanatRetinoMotionLoc_s03xx20150324.m 
%                   by: Steeve Laquitaine
%                 date: 15/03/21
%              Purpose: An fMRI experiment with 3D anatomy, session anatomy and retinotopy 
%            subjectID: s03xx
%                Usage: 
%
%                       run each line for each fMRI scan
%
%note: Requires mgl and mrTools distributions to work.
%
%--------
%Overview
%--------
%We chose MUX2 ARC1 fMRI sequence because we had the most signal. 
%Early visual areas voxels showed the highest coherence values
%(Correlation analysis) for MUX2 ARC1
%
%
% Total time: 40 min
%
% 1) CNI 3D  canonical anatomy .9 mm (1 scan, 4-5 min)
%..................................................
% Canonical anatomy (GE standard) 0.9 mm iso voxel resolution
%
%
%
% 2) 1 CNI 3D session anatomy 1.2 mm (1 scan, 2 min)
%........................................................
%We tested 5 planes to check that the file headers work. They do.
%We will run a coronal plane that covers early visual and LIP
%
% - coronal plane
%
%
% 
% 3) fMRI retinotopy (8 scans, 32 min)
%.....................................
%
% Traveling-wave retinotopic mapping experiment and fixation task
% A scan is a stimulus ring that revolves over the full 
% circle (a cycle) over 24 sec 10 times (10 cycles lasts 4 min) 
%
%We do not run the retinotopy with the projector goggle-like mask because
%we still have things to debug (e.g., black screen should be grey in between run)
%
% e.g. parameters for TR = 1.5s
% - 16 volumes per cycle (= cycle period./TR  [period = 24 sec; TR = 1.5 sec])
% - 168 volumes per scan (10 cycles * 16 vols + 8 vols thrown out + n MUX vols)
% 
%
%How to calculate the scan parameters
%
%
%period
%------
% If you set 'stimulusPeriod=24' and TR = 1.4s
% Actual period = ceil(24s/1.4s)*1.4s) = 25.2s
%
%number of volumes
%-----------------
% If you set 'numCycles=10','initialHalfCycle=1', TR=1.4s and 'stimulusPeriod=24'
% Number of vols = ceil(24s/1.4s) vols/cycles * 10 cycles + .5*ceil(24s/1.4s) vols/cycles  
%
%
%
%
% The code requires .mglScreenParams.mat and saves a .mat file for each scan 
% in a folder ~data/mglRetinotopy_proj with
% - stimulus parameters
% - screen parameters
% - timing of the volumes acquired by the fMRI scanner
%       
%to test simulate a scan run:
%
%            mglSimulateRun(1.57,176,3,'a'); 
%            gruRetinotopy('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
%
%
%
%COMMAMDS
%--------
%
% 1) 1 canonical anatomy : T1w .9mm (CNI software)
%  - save Rx 
%  - run scan
%
% 2) 1 session anatomy : T1w 1.2 mm (CNI software)
% - coronal
% 
%
% 3) fMRI (retinotopy, CNI software and command line)
%
%
%CCW wedges (4:24 min)
%we set: 10 cycles + half cycle , 24s/cycle   , 18 vols/cycle + 9 vols, MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR , 189 vols 
%we get: 10 cycles + half cycle , 25.2s/cycle , 18 vols/cycle + 9 vols, MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR , 189 vols - duration=189*1.4s=264.6s
%mglSimulateRun(1.4,189,1,'a');
gruRetinotopy('displayName=cniLCD','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');

%CW wedges (4:24 min)
gruRetinotopy('displayName=cniLCD','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');

%expanding rings (4:24 min)
gruRetinotopy('displayName=cniLCD','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');

%contracting rings (4:24 min)
gruRetinotopy('displayName=cniLCD','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');

%CCW wedges (4:24 min)
gruRetinotopy('displayName=cniLCD','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');

%CW wedges (4:24 min)
gruRetinotopy('displayName=cniLCD','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');


%-----
%BARS
%-----
%bars (4:49 min)
%we set: 10 cycles + half cycle/blank (3), 24s/cycle   , 18 vols/cycle + 3*9 vols (blanks), MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR ,  207 vols 
%we get: 10 cycles + half cycle/blank (3), 25.2s/cycle , 18 vols/cycle + 3*9 vols (blanks), MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR ,  207 vols , duration=207*1.4=289.8s
mglSimulateRun(1.4,207,1,'a');
gruRetinotopy('displayName=cniLCD','bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','blanks=3','fixedRandom=1','doEyeCalib=-1','initialHalfCycle=1');

gruRetinotopy('displayName=cniLCD','bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','blanks=3','fixedRandom=1','doEyeCalib=-1','initialHalfCycle=1');

gruRetinotopy('displayName=cniLCD','bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','blanks=3','fixedRandom=1','doEyeCalib=-1','initialHalfCycle=1');

gruRetinotopy('displayName=cniLCD','bars=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','blanks=3','fixedRandom=1','doEyeCalib=-1','initialHalfCycle=1');


%17 vols/ cycle - 8.5 (8?) vols half cycle
%mglSimulateRun(1.4,178,10,'a');
%gruRetinotopy('displayName=smallProjector','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');
