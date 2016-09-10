%weijima.m
%
%        $Id: weijima.m
%      usage: weijima
%         by: justin gardner, dylan cable, steeve laquitaine
%       date: 09/07/06, 04/?/15
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%    purpose: weijima
%      usage:
%           
%           weijima('displayName=fMRIprojFlex')
%
%
%Inputs: 
%       'displayName': set display parameters
%                      choose a display available in mglEditScreenParams              
%               e.g., 'displayName=fMRIprojFlex'             
%
%
%Current description:
%         default is to replicate weijima fMRI experiment Nature
%         Neuroscience, 2015
%         
%         An fMRI run is structured: 
%             - start with 4s fixation
%             - 18 x [1.5s grating - 6.5s  fixation - 4s response - 4s fixation]
%             - ends with 4s fixation
%
%       ?parameters
%            - grating (20 deg, 2.5 cycle/deg,counterphasing)
%            - inner 
%            - 4 contrasts(0.012, 0.0145, 0.0165, 0.03)
%            - 72 orientation?(trial-distribution defined by a discrete 
%              vonMises(5:5:360;35,33.33) over ~200 trials (192 exactly))
%
%
%Weijima
%         An fMRI run is structured: 
%             - start with 4s fixation
%             - 18 x [1.5s grating - 6.5s  fixation - 4s response - 4s fixation]
%             - ends with 4s fixation
%
%       ?parameters
%            - grating (20 deg, 2.5 cycle/deg,counterphasing)
%            - inner/outer radius : 1.5/7.5 deg
%            - randomize spatial phase
%            - 2 Hz spatial frequency
%            - linearly decreasing contrast from outer radius out
%              and form inner radius in (becomes 0 at 0.5 degree 
%              radius in and ...out)
%            - 1 contrast (0.1)
%            - 72 orientations?(trial-distribution is discrete uniform
%              5:5:360 over ~200 trials
%
function myscreen = weijima(varargin)

%get args
displayName = 'fMRIprojFlex';
eval(evalargs(varargin,0,0,{'displayName'}));

%parameters
global GRATING_WIDTH GRATING_HEIGHT SF SDX SDY PHASE_CHANGE_PERIOD INITIAL_X INITIAL_Y SCREEN_WIDTH SCREEN_HEIGHT SELECTION_SIZE NUM_TRIALS NUM_SAMPLES STIMLEN REPLEN FEEDLEN MASKLEN RMAX RMIN MIN_ANGLE SCREEN_NUMBER REPORT_INNER_R REPORT_OUTER_R COMPLETE_FEEDBACK CONTRAST MASK_CONTRAST;
GRATING_WIDTH  = 20;
GRATING_HEIGHT = 20;
SELECTION_SIZE = 5;
SF             = 2.5;        %num cycles
SDX            = 1;
SDY            = 1;
PHASE_CHANGE_PERIOD = 0.200; %(counterphasing) change phase every 200 ms
NUM_TRIALS     = 200;
NUM_SAMPLES    = NUM_TRIALS;
STIMLEN        = 2.00; %300 ms
CONTRAST       = [0.0125 0.0145 0.0165 0.0300];
REPLEN         = 5;
MASKLEN        = 1.5;
MASK_CONTRAST  = 0.5;
FEEDLEN        = 0.2;
MIN_ANGLE      = 5;
REPORT_INNER_R = 1.5;
REPORT_OUTER_R = 7.5;
COMPLETE_FEEDBACK = false;
averageAngle   = 35;
concentrationParam = 33.33;
orientations   = 5:5:360;

%orientation trial-distribution defined by a vonMises
myPdf = vmPdfs(orientations,averageAngle,concentrationParam,'norm');
discrete_dist = genDiscreteDist(myPdf,NUM_TRIALS/length(CONTRAST));
displayOrientations = [];
orientationsUsed = [];
for i = 1:length(orientations)
	displayOrientations = [displayOrientations repmat(orientations(i),1,discrete_dist(i))];
	if(discrete_dist(i) > 0)
		orientationsUsed = [orientationsUsed orientations(i)];
	end
end
displayOrientations = displayOrientations(randperm(length(displayOrientations)));
NUM_TRIALS = length(displayOrientations)*length(CONTRAST);

%[tommy Fs] = audioread('eat.mp3');
%sound(tommy,Fs);

% if TR == .75
%     if ~mglGetParam('ignoreInitialVols')==16
%         warning('mux8 script: ignoreInitialVols was set incorrectly, setting to 16');
%         mglSetParam('ignoreInitialVols',16);
%     end
% elseif TR == 1.4
%     if ~mglGetParam('ignoreInitialVols')==4
%         warning('mux2 script: ignoreInitialVols was set incorrectly, setting to 4');
%         mglSetParam('ignoreInitialVols',4);
%     end


% initalize the screen
% myscreen.autoCloseScreen = 1;
% myscreen.displayname = 'projector';
% myscreen.displayname = displayName;
%myscreen.screenNumber=SCREEN_NUMBER;
%myscreen.fixPosX = 0;
%myscreen.fixPosY = 0;
%save data in user dir data/SLtaskLocEst
%myscreen.saveData = 1;
%home = getenv('HOME');
%mkdir([home '/data/dylan'])
%myscreen.datadir = [home '/data/dylan'];
%myscreen.stimfile = [home '/data/dylan'];

%set display parameters
myscreen.displayname = displayName;

%open screen
myscreen = initScreen(myscreen);

%hide mouse
mglDisplayCursor(0);

%clear screen
mglClearScreen(0.5);
mglFlush;
mglClearScreen;

SCREEN_WIDTH = mglGetParam('screenWidth');
SCREEN_HEIGHT = mglGetParam('screenHeight');
INITIAL_X = SCREEN_WIDTH/2 ;
INITIAL_Y = SCREEN_HEIGHT/2;
disp(SCREEN_WIDTH);
disp(SCREEN_HEIGHT);
screenSize = min(SCREEN_WIDTH,SCREEN_HEIGHT);
%RMAX = screenSize/2 * 0.75;
%RMIN = screenSize/2 * 0.1;

%initialize grating orientations
%initAngleSequence = round(rand(NUM_TRIALS,1)*360+0.5);

%a top-up period of the same direction
task{1}{1}.segmin = [0.2 STIMLEN MASKLEN REPLEN FEEDLEN];
task{1}{1}.segmax = [0.2 STIMLEN MASKLEN REPLEN FEEDLEN];
task{1}{1}.getResponse = [0 0 0 1 0];
task{1}{1}.parameter.contrast = CONTRAST;
task{1}{1}.parameter.dir = displayOrientations;
task{1}{1}.random = 1;
task{1}{1}.numTrials = NUM_TRIALS;


% task{1}{1}.waitForBacktick = 1;

%initialize variables to analyse
task{1}{1}.randVars.uniform.initAngleDeg = [1:360];
task{1}{1}.randVars.calculated.meanResponseRecorded = nan;
task{1}{1}.randVars.calculated.meanResponse = nan;
task{1}{1}.randVars.calculated.meanResponseTime = nan;
task{1}{1}.randVars.calculated.phaseSequence = {nan};
task{1}{1}.randVars.calculated.angleReport = nan;
task{1}{1}.randVars.calculated.confidence = nan;
task{1}{1}.randVars.calculated.angleStart = nan;
task{1}{1}.randVars.calculated.angleLen = nan;
task{1}{1}.randVars.calculated.correct = nan;

%store parameters
task{1}{1}.data.concentrationParam = concentrationParam;
task{1}{1}.data.averageAngle = averageAngle;
task{1}{1}.data.orientationsUsed = orientationsUsed;
task{1}{1}.data.REPORT_INNER_R = REPORT_INNER_R;
task{1}{1}.data.REPORT_OUTER_R = REPORT_OUTER_R;
task{1}{1}.data.COMPLETE_FEEDBACK = COMPLETE_FEEDBACK;

% initialize our task
for phaseNum = 1:length(task{1})
[task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);
end

% init the stimulus
clear global stimulus;
global stimulus Inicoord INITIAL_X INITIAL_Y RMAX RMIN;
myscreen = initStimulus('stimulus',myscreen);
%myscreen.datadir = '~/Documents/Stimulus';

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%% Colors
stimulus.colors.rmed = 127.75;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.

stimulus.colors.reservedBottom = [0 0 0; 1 1 1; repmat(0.5 - MASK_CONTRAST/2,1,3)]; % fixation cross colors
stimulus.colors.reservedTop = [repmat(0.5 + MASK_CONTRAST/2,1,3); 1 0 0; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = [0/255 0/255 0/255]; stimulus.colors.white = [1/255 1/255 1/255]; stimulus.colors.darkGray = 2/255;
stimulus.colors.red = [254/255 254/255 254/255]; stimulus.colors.green = [255/255 255/255 255/255]; stimulus.colors.lightGray = 253/255;
stimulus.colors.nReserved = 3; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Gamma Table Initialization
stimulus = initPatch(stimulus,myscreen);
stimulus.type = 'patch';
Inicoord.x = INITIAL_X;
Inicoord.y = INITIAL_Y;
Inicoord.rmax = RMAX;
Inicoord.rmin = RMIN;

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
% update the dots
[task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
% flip screen
myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus STIMLEN PHASE_CHANGE_PERIOD INITIAL_X INITIAL_Y SCREEN_NUMBER;
mglSetMousePosition(INITIAL_X,INITIAL_Y,SCREEN_NUMBER);
if(task.thistrial.thisseg == 1)
%initialize task variables
task.thistrial.phaseSequences = zeros(1,round(1+STIMLEN/PHASE_CHANGE_PERIOD));
task.thistrial.phaseIndex = 1;
end


if (task.thistrial.thisseg == 2)
stimulus.patch.contrast = task.thistrial.contrast;
stimulus.patch.angle = task.thistrial.dir;
stimulus.patch.phaseChangeTime = 0;
stimulus.patch.startTime = task.thistrial.segstart;
end
if(task.thistrial.thisseg == 3)
	stimulus.mask.isSet = 0;
	stimulus.mask.tex = 0;
end
if (task.thistrial.thisseg == 4)
task.thistrial.phaseSequence = stimulus.phaseSequence;
task.thistrial.meanResponseRecorded = 0;
end
if (task.thistrial.thisseg == 5)
task.thistrial.confidenceResponseRecorded = 0;
end
stimulus.patch.dir = task.thistrial.dir;
stimulus.response.initangle = task.thistrial.initAngleDeg;
%stimulus = updatePatch(stimulus,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus MIN_ANGLE REPORT_INNER_R REPORT_OUTER_R COMPLETE_FEEDBACK;
mglClearScreen;

%motion
if(task.thistrial.thisseg == 2)
	setGammaTable_flowMax(task.thistrial.contrast);
	stimulus = updatePatch(stimulus,myscreen);
end

if(task.thistrial.thisseg == 3)
	%setGammaTable_flowMax(1);
	stimulus = updateMask(stimulus,myscreen);
end
%confidence response
if (task.thistrial.thisseg == 4)

%draw resp arrow only foveal (separated from peripheral motion)


%update until response
%task.thistrial.gotResponse is 1 for keyboard
%task.thistrial.gotResponse is 2 for mouse click
if task.thistrial.meanResponseRecorded == 0
%Position mouse on aperture
%mi = mglGetMouse(myscreen.screenNumber); %(pixels)
mi = mglGetMouse();
mouseinfo.x.pix = mi.x;
mouseinfo.y.pix = mi.y;

%Position mouse on aperture by moving mouse horizontally
mouseinfo = Pos2apByMovingHoriz(mouseinfo,stimulus,task);
%mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo);

%Position mouse on aperture
%mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo);

%Interface PowerMate
%mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);

%draw a white arrow (radius)
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.black,60);
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.black,60);
%mglGluDisk(-mouseinfo.x.angle,-mouseinfo.y.angle,0.2,0);
task.thistrial.angleReport = 90 - mouseinfo.angle; %wierd coordinates
task.thistrial.confidence = mouseinfo.conf;
task.thistrial.angleStart = 90 - mouseinfo.angleStart;
task.thistrial.angleLen = mouseinfo.angleLen;
%mglLines2(100,100,mouseinfo.x.angle,mouseinfo.y.angle,1,[1 1 1],1);

else
end
end
if(task.thistrial.thisseg == 5)
if task.thistrial.meanResponseRecorded == 1
mouseinfo.angleLen = task.thistrial.angleLen;
mouseinfo.angleStart = 90 - task.thistrial.angleStart; % for drawring
myDir = 90 - task.thistrial.dir;

%draw a white arrow (radius)
%EXPLICIT FEEDBACK
%if(mod(task.thistrial.dir - (task.thistrial.angleStart - MIN_ANGLE/2),180) < task.thistrial.angleLen + MIN_ANGLE)
%mglGluPartialDisk(0,0,2.9,3,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2,mouseinfo.angleLen + MIN_ANGLE,[0 1 0],60);
%mglGluPartialDisk(0,0,2.9,3,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2 + 180,mouseinfo.angleLen + MIN_ANGLE,[0 1 0],60);
%task.thistrial.correct = 1;
%else
%mglGluPartialDisk(0,0,2.9,3,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2,mouseinfo.angleLen + MIN_ANGLE,[1 0 0],60);
%mglGluPartialDisk(0,0,2.9,3,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2 + 180,mouseinfo.angleLen + MIN_ANGLE,[1 0 0],60);
%task.thistrial.correct = 0;
%end
%draw a white arrow (radius)
if(mod(myDir + MIN_ANGLE/2 - mouseinfo.angleStart,180) < task.thistrial.angleLen + MIN_ANGLE)
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.green,60);
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.green,60);
if(COMPLETE_FEEDBACK)
	mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.white,60);
	mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.white,60);
end
task.thistrial.correct = 1;
else
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.red,60);
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.red,60);
if(COMPLETE_FEEDBACK)
	mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.white,60);
	mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,myDir - (MIN_ANGLE + mouseinfo.angleLen)/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.white,60);
end
task.thistrial.correct = 0;
end


else
mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,360,stimulus.colors.red,60);
task.thistrial.correct = 0;
end
end
if(task.thistrial.thisseg <= 5)
mglGluDisk(0,0,0.2,0);      %fixation (deg)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the orientation stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initPatch(stimulus,myscreen)
global GRATING_WIDTH GRATING_HEIGHT SF SDX SDY PHASE_CHANGE_PERIOD SELECTION_SIZE STIMLEN;
% convert the passed in parameters to real units
if ~isfield(stimulus,'patch') || ~isfield(stimulus.patch,'width'), stimulus.patch.width = GRATING_WIDTH;,end
if ~isfield(stimulus.patch,'height'), stimulus.patch.height = GRATING_HEIGHT;,end
if ~isfield(stimulus.patch,'sf'), stimulus.patch.sf = SF;,end
if ~isfield(stimulus.patch,'sdx'), stimulus.patch.sdx = SDX;,end
if ~isfield(stimulus.patch,'sdy'), stimulus.patch.sdy = SDY;,end
if ~isfield(stimulus.patch,'angle'), stimulus.patch.angle = 0;,end
if ~isfield(stimulus.patch,'phaseChangeTime'), stimulus.patch.phaseChangeTime = 0;,end %timeToChangePhase
if ~isfield(stimulus.patch,'phaseChangePeriod'), stimulus.patch.phaseChangePeriod = PHASE_CHANGE_PERIOD;,end %timeToChangePhase
if ~isfield(stimulus.patch,'startTime'), stimulus.patch.startTime = 0;,end
if ~isfield(stimulus.patch,'texture'), stimulus.patch.texture = [];,end
if ~isfield(stimulus,'response') || ~isfield(stimulus.response,'rmax'), stimulus.response.rmax = SELECTION_SIZE;,end
if ~isfield(stimulus.response,'initangle'), stimulus.response.angle = 0;,end
if ~isfield(stimulus,'trialIndex'), stimulus.trialIndex = 0;,end
stimulus.phaseindex = 1;
stimulus.phaseSequence = cell(1,round(1+STIMLEN/PHASE_CHANGE_PERIOD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to update patch poisition and draw to screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = updatePatch(stimulus,myscreen)
if(stimulus.patch.contrast > 1E-6)
	if(mglGetSecs(stimulus.patch.startTime) > stimulus.patch.phaseChangeTime)
		if(~ (isempty(stimulus.patch.texture)))
			mglDeleteTexture(stimulus.patch.texture);
		end
		phase = rand()*360;
		if(stimulus.phaseindex <= length(stimulus.phaseSequence))
			stimulus.phaseSequence{stimulus.phaseindex} = phase;
			stimulus.phaseindex = stimulus.phaseindex + 1;
		end
		grating = mglMakeGrating(stimulus.patch.width,stimulus.patch.height,stimulus.patch.sf,stimulus.patch.angle,phase);
		gaussian = mglMakeGaussian(stimulus.patch.width,stimulus.patch.height,stimulus.patch.sdx,stimulus.patch.sdy);
		gabor = convertContrast(255*(grating.*gaussian+1)/2);
		tex = mglCreateTexture(gabor);
		mglBltTexture(tex,[0 0]);
		stimulus.patch.texture = tex;
		stimulus.patch.phaseChangeTime = mglGetSecs(stimulus.patch.startTime) + stimulus.patch.phaseChangePeriod;
		else
		mglBltTexture(stimulus.patch.texture,[0 0]);
	end
end

function stimulus = updateMask(stimulus,myscreen)
if(~ stimulus.mask.isSet)
	stimulus.mask.isSet = 1;
	sigma = 40; %parameter for gaussian
	%make white noise texture
	gaussian = mglMakeGaussian(stimulus.patch.width,stimulus.patch.height,stimulus.patch.sdx,stimulus.patch.sdy);
	L = length(gaussian)/2;
	startInd = round(L/2);
	endInd = round(L + L/2);
	gaussian = gaussian(startInd:endInd,startInd:endInd);
	grating = zeros(size(gaussian,1),size(gaussian,2));
	for i = 1:size(grating,1)
		for j = 1:size(grating,2)
			grating(i,j) =  rand(); %norminv(rand(),0,1); %min(max(normSampleScalar(255/2,50),0),255);
		end
	end
	grating = mat2gray(grating);
	freq = fftshift(fft2(grating));
	for i = 1:size(grating,1)
	    for j = 1:size(grating,2)
	        d = sqrt((i - size(grating,1)/2)^2 + (j - size(grating,2)/2)^2);
	        weight = 1/sqrt(2*sigma)*exp(-d^2/(2*sigma^2));
	        freq(i,j) =  weight*freq(i,j); %min(max(normSampleScalar(255/2,50),0),255);
	    end
	end
	freq = fftshift(freq);
	grating = abs(ifft2(freq));
	Y = 200*mat2gray(grating);
	average = mean(mean(Y));
	Y = min(max(Y  - average,-255/2),255/2).*gaussian + 255/2;
	threshold1 = 127.0;
	threshold2 = 128.0;
	%Splits the image into either lightGray, medium, or darkGray
	A = (Y > threshold2) *stimulus.colors.darkGray*255;
	B = (Y < threshold1) *stimulus.colors.lightGray*255;
	C = ((Y > threshold1) & (Y < threshold2))*127;
	D = A + B + C;
	for i = 1:size(Y,1)
		break
		for j = 1:size(Y,2)
			if(Y(i,j) > threshold1)
				if(Y(i,j) > threshold2)
					Y(i,j) = stimulus.colors.darkGray*255;
				else
					Y(i,j) = 127;
				end
			else
				Y(i,j) = stimulus.colors.lightGray*255;
			end
		end
	end
	tex = mglCreateTexture(D);
	mglBltTexture(tex,[0 0]);
	stimulus.mask.tex = tex;
else
	mglBltTexture(stimulus.mask.tex,[0 0]);
end
	

%function to convert from PowerMate's to aperture's coordinates
function mouseinfo = Pos2apByMovingHoriz(mouseinfo,stimulus,task)
global Inicoord maxX SCREEN_WIDTH;
if(task.thistrial.thisseg == 4)
	scrLen = SCREEN_WIDTH;
	x = mouseinfo.x.pix - Inicoord.x - scrLen;
	if(x > maxX)
		maxX = x
	end
	conf = .95;
	arclen = (1 - conf) * 180; % in degrees
	mouseinfo.angle = (x/scrLen + 0.5)*360 + stimulus.response.initangle;
	mouseinfo.angleStart = mouseinfo.angle - arclen/2;
	mouseinfo.angleLen = arclen;
	mouseinfo.conf = conf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called when subject responds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task, myscreen)

% get correct or incorrect
response = find(task.thistrial.buttonState);
if(task.thistrial.thisseg == 4 && task.thistrial.meanResponseRecorded == 0)
task.thistrial.meanResponseRecorded = 1;
task.thistrial.meanResponse = response;
task.thistrial.meanResponseTime = mglGetSecs(task.thistrial.segStartSeconds);
task = jumpSegment(task);
end

%returns between 0 and 180
function sample = normSample(average, variance)

N = 1000; %num_sum
x = sum(rand(1,N) - 0.5)/sqrt(N/12); %0/1 random
sample = mod(x*sqrt(variance) + average,180);

%returns a scalar, not an angle
function sample = normSampleScalar(average, variance)

N = 1000; %num_sum
x = sum(rand(1,N) - 0.5)/sqrt(N/12); %0/1 random
sample = x*sqrt(variance) + average;

function setGammaTable_flowMax(maxContrast)

global stimulus;

%Correct nan values in gamma table
if(isnan(stimulus.linearizedGammaTable.redTable(1)))
	stimulus.linearizedGammaTable.redTable(1) = 0;
end
if(isnan(stimulus.linearizedGammaTable.greenTable(1)))
	stimulus.linearizedGammaTable.greenTable(1) = 0;
end
if(isnan(stimulus.linearizedGammaTable.blueTable(1)))
	stimulus.linearizedGammaTable.blueTable(1) = 0;
end
% set the bottom. Note: assumes all three linearizedGammaTables are identical.
gammaTable = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,stimulus.colors.reservedBottom,'linear');

% set the gamma table
if maxContrast == 1
    % create the rest of the gamma table
    cmax = 1;cmin = 0;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
elseif maxContrast > 0
    % create the rest of the gamma table
    cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nUnreserved-1)):cmax;

    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
else
    % if we are asked for 0 contrast then simply set all the values to gray
    redLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    greenLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
    blueLinearized = repmat(.5,1,stimulus.colors.nUnreserved);
end

% add to the table!
gammaTable((stimulus.colors.mrmin:stimulus.colors.mrmax)+1,:)=[redLinearized;greenLinearized;blueLinearized]';

% set the top
gammaTable = [gammaTable; interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,stimulus.colors.reservedTop,'linear')];
if size(gammaTable,1)~=256
    disp('(setGammaTable) Failure: Incorrect number of colors in gamma table produced');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.curMaxContrast = maxContrast;

%converts contrast to gamma table index according to 
function contrast = convertContrast(contrast)
global stimulus;
contrast = stimulus.colors.nReserved + contrast*(stimulus.colors.nUnreserved*1.0/256);



