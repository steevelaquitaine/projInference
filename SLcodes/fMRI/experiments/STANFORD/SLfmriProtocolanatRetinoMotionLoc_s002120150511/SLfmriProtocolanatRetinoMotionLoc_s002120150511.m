
% SLfmriProtocolanatRetinoMotionLoc_s002120150511.m
% 
%                  $Id: SLfmriProtocolanatRetinoMotionLoc_s002120150511.m 
%                   by: Steeve Laquitaine
%                 date: 15/05/11
%              Purpose: An fMRI experiment with 3D anatomy, session anatomy and retinotopy 
%            subjectID: s0021
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
% - 48 volumes per cycle (= cycle period./TR  [period = 24 sec; TR = 1.4 sec])
% - 488 volumes per scan (10 cycles * 48 vols + 8 vols junkFrames + n MUX vols)
%
%Also, calibration TRs are always acquired by the CNI scanner to help image
%reconstruction. The number of calibration TRs depends on your sequence (Mux and Arc)
%   total number of TRs acquired = n + calibration TRs
%   number of calibration TRs = mux * arc * num_mux_cycle (CNI wiki, num_mux_cycle = 2 by default)
%   number of TRs in nifti = n + num_mux_cycle
%   number in the phases per location field: n + mux * n_mux_cycle
%
%   total number of TRs acquired = 480 + calibration TRs
%   number of calibration TRs = 2 * 1 * 2 = 4
%   number of TRs in nifti = 480 + 2 = 482
%   number in the phases per location field: 480 + 2 * 2 = 484 (what you enter in soft)
%
%!! Always make sure that you have 10 cycles of usable TRs (480 TRs) in the
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
% If you set 'numCycles=10','initialHalfCycle=0', TR=1.4s and 'stimulusPeriod=24'
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
%            mglSimulateRun(1.57,486,3,'a'); 
%            gruRetinotopy('displayName=smallProjector','rings=1','direction=1','numCycles=10','stimulusPeriod=24','synchToVolEachCycle=0','doEyeCalib=0','initialHalfCycle=0');
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
%for mux 2 ARC1 we set: 10 cycles + half cycle , 24s/cycle , 48 vols/cycle + 8 vols (2 are thrown out by mux 2, 6 are junkframes), 2.4 mm iso vox.res, 1.4s TR , 488 vols + mux2*numux = 182 vols
%   10*48+floor(48/2) + 4 = 182 vols
%
%for mux 8 ARC1 we set: 10 cycles, 24s/cycle   , 48 vols/cycle, 2.4 mm iso vox.res, 0.5s TR, 480 vols + mux8*nummux = 496 vols
%   10*48 + 16 = 496 vols

% mglPostEvent('quit')
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%CW wedges (182 vols alias "phases per location", 4:14 min)
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%expanding rings (182 vols alias "phases per location", 4:14 min)
mglRetinotopy('displayName=fMRIproj16','rings=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%contracting rings (182 vols alias "phases per location", 4:14 min)
mglRetinotopy('displayName=fMRIproj16','rings=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%CCW wedges (182 vols alias "phases per location", 4:14 min)
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');

%CW wedges (182 vols alias "phases per location", 4:14 min)
mglRetinotopy('displayName=fMRIproj16','wedges=1','direction=-1','numCycles=10','stimulusPeriod=24','stepsPerCycle',48,'doEyeCalib=0','initialHalfCycle=0');


%-----
%BARS
%-----
%for mux 2 arc 1 we set: 8 directions of 48 steps (a cycle) + 3 blanks of 1/2 cycle (8 steps), period: 24s/cycle , ~48 vols/cycle + 3*8 vols(blanks), 2.4 mm iso vox.res , 1.4s TR , 160 vols + mux2*nummux = 164 vols (3:49 min)
%for mux 8 arc 1 we set: same with 0.5s TR , 48 vols/cycle , 8*48 + 3*24 = 456 vols + mux8*nummux = 472 vols (3:55 min)


%What the code does: a cycle is when the bars goes from one side of the screen 
%to the other in the same direction. The cycle consists of 48 steps and last 
%24s which is ~ 48 vols for a Tr = 1.4s (24s/1.4s). The bars always go through 
%8 directions ([0:45:359]) and "numCycles" tag is useless here. You should
%add "blanks" which are period of grey screen (no stim).
%Total number of vols will be 
%
%           numVols = 8 dir * floor(24s/1.4s) + 3 blanks * floor(floor(24s/1.4s)/2) half cycles blanks
%
%The bars code essentially requires "stimulusPeriod" and "stepsPerCycle" or
%"synchToVolEachCycle".

%1 (164 vols alias "phases per location")
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%2 (164 vols alias "phases per location")
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%3 (164 vols alias "phases per location")
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');

%4 (164 vols alias "phases per location")
mglRetinotopy('displayName=fMRIproj16','bars=1','fixedRandom=1','stimulusPeriod=24','stepsPerCycle',48,'blanks=3','doEyeCalib=0');












