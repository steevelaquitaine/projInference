
% SLfmriProtocolanatRetino_s030020150403.m
% 
%                  $Id: SLfmriProtocolanatRetino_s030020150403.m 
%                   by: Steeve Laquitaine
%                 date: 15/04/03
%              Purpose: An fMRI experiment with 3D anatomy, session anatomy and retinotopy 
%            subjectID: s0300
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
% - 17 volumes per cycle (= cycle period./TR  [period = 24 sec; TR = 1.4 sec])
% - 178 volumes per scan (10 cycles * 17 vols + 8 vols junkFrames + n MUX vols)
%
%Also, calibration TRs are always acquired by the CNI scanner to help image
%reconstruction. The number of calibration TRs depends on your sequence (Mux and Arc)
%   total number of TRs acquired = n + calibration TRs
%   number of calibration TRs = mux * arc * num_mux_cycle (CNI wiki, num_mux_cycle = 2 by default)
%   number of TRs in nifti = n + num_mux_cycle
%   number in the phases per location field: n + mux * n_mux_cycle
%
%   total number of TRs acquired = 170 + calibration TRs
%   number of calibration TRs = 2 * 1 * 2 = 4
%   number of TRs in nifti = 170 + 2 = 172
%   number in the phases per location field: 170 + 2 * 2 = 174 (what you enter in soft)
%
%!! Always make sure that you have 10 cycles of usable TRs (170 TRs) in the
%end for the wedges and rings correlation analysis.!!
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
%to test simulate a scan run in debug mode :
%
%            mglSimulateRun(1.57,176,3,'a'); 
%            gruRetinotopy('displayName=smallProjector','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=1','doEyeCalib=-1','initialHalfCycle=1');
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
% 3) fMRI (retinotopy, CNI softw. and matlab)
%
%
%CCW wedges (4:9 min)
%we set: 10 cycles + half cycle , 24s/cycle , 17 vols/cycle + 8 vols (2 are thrown out by mux 2, 6 are junkframes), MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR , 178 vols 
gruRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%CW wedges (178 vols alias "phases per location", 4:9 min)
gruRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%expanding rings (178 vols alias "phases per location", 4:9 min)
gruRetinotopy('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%contracting rings (178 vols alias "phases per location", 4:9 min)
gruRetinotopy('displayName=fMRIproj16','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%CCW wedges (178 vols alias "phases per location", 4:9 min)
gruRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');

%CW wedges (178 vols alias "phases per location", 4:9 min)
gruRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',17,'doEyeCalib=-1','initialHalfCycle=1');


%-----
%BARS
%-----
%we set 8 directions of 17 steps (a cycle) + 3 blanks of 1/2 cycle (8 steps), period: 24s/cycle , ~17 vols/cycle + 3*8 vols(blanks), MUX2-ARC1 , 2.4 mm iso vox.res , 1.4s TR , 160 vols, duration: 224 s (3:43s)

%What the code does: a cycle is when the bars goes from one side of the screen 
%to the other in the same direction. The cycle consists of 17 steps and last 
%24s which is ~ 17 vols for a Tr = 1.4s (24s/1.4s). The bars always go through 
%8 directions ([0:45:359]) and "numCycles" tag is useless here. You should
%add "blanks" which are period of grey screen (no stim).
%Total number of vols will be 
%
%           numVols = 8 dir * floor(24s/1.4s) + 3 blanks * floor(floor(24s/1.4s)/2) half cycles blanks
%
%The bars code essentially requires "stimulusPeriod" and "stepsPerCycle" or
%"synchToVolEachCycle".

%1 (160 vols alias "phases per location")
%mglSimulateRun(1.4,8*17+3*floor(17/2),10,'a');
gruRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');

%2 (160 vols alias "phases per location")
gruRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');

%3 (160 vols alias "phases per location")
gruRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');

%4 (160 vols alias "phases per location")
gruRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',17,'blanks=3','doEyeCalib=-1');












