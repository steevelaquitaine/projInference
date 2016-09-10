%weijima2.m
%
%        $Id: weijima.m
%      usage: weijima
%         by: justin gardner, dylan cable, steeve laquitaine
%       date: 09/07/06, 04/?/15
%  copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%     status: complete
%    purpose: Replicate weijima fMRI experiment Nature
%             Neuroscience, 2015
%
%      usage:
%           
%           mglSimulateRun(0.5,608,1,'a') %to simulate scanner
%           weijima2('displayName=fMRIprojFlex')
% 
%
%Inputs:
%
%       'displayName': set display parameters
%                      choose a display aavailable in mglEditScreenParams
%                      e.g., 'displayName=fMRIprojFlex'
%
%outputs:
%
%       %a stimfile saved in ~/data/ that you can load
%       %then retrieve task parameters
%       a = getTaskParameter(myscreen,task)
%       a
%       %get angles reported by subjects
%       a{1}(2).randVars.angleReport  
%
%
%Description:
%
%         Replicate weijima fMRI experiment Nature
%         Neuroscience, 2015
%
%         An fMRI run is structured:
%           - 4sec fixation - 18 x [4s fixation - 1.5s grating - 6.5s
%           ISI - 4s response] - 4 sec fixation
%           - a run lasts 296 sec (4 + 18*(4+1.5+6.5+4) + 4; 4:55min)
%           in stanford cni scanner a run is
%           - 608 vols of .5s TR + (592 task vols + 16 initial vols of mux8)
%           - 304 sec (5:03min)
%           then screen remains gray (press esc)
%
%       stimulus parametersa
%            - grating (spatial frequency: 2.5 cycle/deg)
%            - inner/outer radius : 1.5/7.5 deg
%            - randomize spatialaa phaaaaaaaase (counterphasing)
%            - 2 Hz temporal frequency (phase changes every 0.5s)
%            - linearly decreasing contrast from outer radius out
%              and from inner radius in (becomes 0 at 0.5 degrees
%              radius in and ...out)
%            - 1 contrast (0.1)
%            - 6 orientatioans (trial-distribution is discrete uniform
%              10:30:160 over ~200 trials: different from wjm who sampled
%              randomly because we want repeats)
%
%       Response:
%
%            - subject moves the mouse horizontally to rotate clockwise/
%              counterclockwise (no confirmation press)

function myscreen = weijima2(varargin)

%get args
displayName = 'fMRIprojFlex';
eval(evalargs(varargin,0,0,{'displayName'}));

%global
global GRATING_WIDTH GRATING_HEIGHT PHASE_CHANGE_PERIOD INITIAL_X INITIAL_Y SCREEN_WIDTH SCREEN_HEIGHT SELECTION_SIZE NUM_TRIALS NUM_SAMPLES STIMLEN REPLEN FEEDLEN MASKLEN RMAX RMIN MIN_ANGLE REPORT_INNER_R REPORT_OUTER_R COMPLETE_FEEDBACK CONTRAST MASK_CONTRAST;

%full screen square grating width and height (deg)
%then clipped with an annular masked at outer
%and inner edges
GRATING_WIDTH  = 20;
GRATING_HEIGHT = 20;
NUM_TRIALS     = 18; %18 (WJM,2015)
NUM_SAMPLES    = NUM_TRIALS;
STIMLEN        = 1.5; %1.5s (WJM,2015)
CONTRAST       = 0.1; %10% (WJM,2015)
REPLEN         = 4;
MASKLEN        = 0;
MASK_CONTRAST  = 0.1;
FEEDLEN        = 0.2;
MIN_ANGLE      = 5;
COMPLETE_FEEDBACK = false;
% averageAngle   = 35;
concentrationParam = 33.33;
orientations   = 10:30:160;

%--------------------------- stimulus parameter ------------------------
global stimulus Inicoord
stimulus.patch.inR  = 1.5;   %annulus inner radius (deg, wjm)
stimulus.patch.outR = 7.5;   %annulus outer radius (deg, wjm)
stimulus.patch.sf   = 1;     %grating spatial frequency (cycle/deg, wjm)
PHASE_CHANGE_PERIOD = 0.500; %phases change every 0.5s (temporal frequency: 2Hz, wjm)


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

myscreen.displayName = displayName;
myscreen = initScreen(myscreen);
mglDisplayCursor(0);
mglClearScreen(0.5);
mglFlush;
mglClearScreen;
SCREEN_WIDTH  = mglGetParam('screenWidth');
SCREEN_HEIGHT = mglGetParam('screenHeight');
INITIAL_X     = SCREEN_WIDTH/2 ;
INITIAL_Y     = SCREEN_HEIGHT/2;
disp(SCREEN_WIDTH);
disp(SCREEN_HEIGHT);
%screenSize = min(SCREEN_WIDTH,SCREEN_HEIGHT);
%initial fixation phase
task{1}{1}.seglen    = 4;
task{1}{1}.numTrials = 1;

%a top-up period of the same direction
%add a segment of 0.1 sec just to record
%subject response
task{1}{2}.segmin = [4 STIMLEN 6.5 REPLEN-0.01 0.01];
task{1}{2}.segmax = [4 STIMLEN 6.5 REPLEN-0.01 0.01];
task{1}{2}.synchToVol = [0 0 0 0 1];%fmri pulses syncs at the end of the segment
task{1}{1}.fudgeLastVolume = 1;
task{1}{1}.waitForBacktick = 1;   %fmri pulse sync stim
task{1}{2}.getResponse = [0 0 0 0 1];
task{1}{2}.parameter.contrast = CONTRAST;
task{1}{2}.parameter.dir = orientations;
task{1}{2}.random = 1;
task{1}{2}.numTrials = NUM_TRIALS;
task{1}{2}.randVars.uniform.initAngleDeg = [1:360];
task{1}{2}.randVars.calculated.meanResponseRecorded = nan;
task{1}{2}.randVars.calculated.meanResponse = nan;
task{1}{2}.randVars.calculated.meanResponseTime = nan;
task{1}{2}.randVars.calculated.phaseSequence = {nan};
task{1}{2}.randVars.calculated.angleReport = nan;
task{1}{2}.randVars.calculated.confidence = nan;
task{1}{2}.randVars.calculated.angleStart = nan;
task{1}{2}.randVars.calculated.angleLen = nan;
task{1}{2}.randVars.calculated.correct = nan;
task{1}{2}.data.concentrationParam = concentrationParam;
% task{1}{2}.data.averageAngle = averageAngle;
task{1}{2}.data.REPORT_INNER_R = REPORT_INNER_R;
task{1}{2}.data.REPORT_OUTER_R = REPORT_OUTER_R;
task{1}{2}.data.COMPLETE_FEEDBACK = COMPLETE_FEEDBACK;

%final fixation phase
task{1}{3}.seglen = 4;
task{1}{3}.numTrials = 1;

%initialize our task
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@responseCallback);
end

%---------------------------- init the stimulus ---------------------------
myscreen = initStimulus('stimulus',myscreen);
%myscreen.datadir = '~/Documents/Stimulus';

%get gamma table
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
global stimulus STIMLEN PHASE_CHANGE_PERIOD INITIAL_X INITIAL_Y;
%toc;
%tic;
if myscreen.screenNumber == 0
    myscreen.screenNumber =1;
end
mglSetMousePosition(INITIAL_X,INITIAL_Y,myscreen.screenNumber);

%fixation and grating
if(task.thistrial.thisphase == 2 && task.thistrial.thisseg == 1)
    %initialize task variables
    task.thistrial.phaseSequences = zeros(1,round(1+STIMLEN/PHASE_CHANGE_PERIOD));
    task.thistrial.phaseIndex = 1;
    setGammaTable_flowMax(task.thistrial.contrast);
end

%grating
if (task.thistrial.thisseg == 2)
    stimulus.patch.contrast = task.thistrial.contrast;
    stimulus.patch.angle = task.thistrial.dir;
    stimulus.patch.phaseChangeTime = 0;
    stimulus.patch.startTime = task.thistrial.segstart;
end

%ISI
%might want to set that back to 1 later
if(task.thistrial.thisseg == 3)
    %setGammaTable_flowMax(1); %full contrast
    setGammaTable_flowMax(task.thistrial.contrast); %full contrast
end

%Response
if (task.thistrial.thisseg == 4)
    task.thistrial.phaseSequence = stimulus.phaseSequence;
    task.thistrial.meanResponseRecorded = 0;
    setGammaTable_flowMax(1); %full contrast
    tic; %start counting time
end

%print out response 
if (task.thistrial.thisseg == 5)   

    %print trial number with info and orientation choosen.
    fprintf('%s %s  %s \n','trial','true angle','response')
    fprintf('    %i    %i   %1.03f  \n \n',task.trialnum,...
        task.thistrial.dir,...
        task.thistrial.angleReport)
end



if(task.thistrial.thisphase == 2)
    stimulus.patch.dir = task.thistrial.dir;
    stimulus.response.initangle = task.thistrial.initAngleDeg;
end
%stimulus = updatePatch(stimulus,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = updateScreenCallback(task, myscreen)

global stimulus MIN_ANGLE REPORT_INNER_R REPORT_OUTER_R COMPLETE_FEEDBACK REPLEN;
mglClearScreen;

%motion
if(task.thistrial.thisphase == 2 && task.thistrial.thisseg == 1)
    
end
if(task.thistrial.thisseg == 2)
    stimulus = updatePatch(stimulus,myscreen);
    %mglGluDisk(0,0,1.5,0)  %foveal black disk
end
if(task.thistrial.thisseg == 3)
    
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
        line_length = 2.8;
        radius = line_length/2;
        TIME_TO_DECAY = 1;
        time_left = min(TIME_TO_DECAY,max(REPLEN - toc,0));
        time_weight = 1 - time_left;
        line_width = 0.07; %degrees;
        line_width_pix = (mglGetParam('xDeviceToPixels')+mglGetParam('yDeviceToPixels'))*0.5*line_width; %convert to pixels
        mglLines2(radius*cosd(90-mouseinfo.angle), radius*sind(90-mouseinfo.angle), -radius*cosd(90-mouseinfo.angle), -radius*sind(90-mouseinfo.angle), line_width_pix, convertContrast([0,0,0]+0.5*time_weight*[255,255,255]));
        %mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.black,60);
        %mglGluPartialDisk(0,0,REPORT_INNER_R,REPORT_OUTER_R,mouseinfo.angleStart - MIN_ANGLE/2 + 180,mouseinfo.angleLen + MIN_ANGLE,stimulus.colors.black,60);
        %mglGluDisk(-mouseinfo.x.angle,-mouseinfo.y.angle,0.2,0);
        task.thistrial.angleReport = 90 - mouseinfo.angle;        %wierd coordinates        
        
        %angle in deg reported by subject
        %always between 1 and 180 deg
        task.thistrial.angleReport = slKeepAngleBetween0to179(task.thistrial.angleReport);
        
        %confidence
        task.thistrial.confidence = mouseinfo.conf;
        task.thistrial.angleStart = 90 - mouseinfo.angleStart;
        task.thistrial.angleLen = mouseinfo.angleLen;
        %mglLines2(100,100,mouseinfo.x.angle,mouseinfo.y.angle,1,[1 1 1],1);
    end    
end
if(true)
    mglGluDisk(0,0,0.25,0,24);      %fixation (deg)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the orientation stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initPatch(stimulus,myscreen)
global GRATING_WIDTH GRATING_HEIGHT SF PHASE_CHANGE_PERIOD SELECTION_SIZE STIMLEN;
% convert the passed in parameters to real units
if ~isfield(stimulus,'patch') || ~isfield(stimulus.patch,'width'), stimulus.patch.width = GRATING_WIDTH;,end
if ~isfield(stimulus.patch,'height'), stimulus.patch.height = GRATING_HEIGHT;,end
if ~isfield(stimulus.patch,'inR'), stimulus.patch.height = 1.5;,end
if ~isfield(stimulus.patch,'outR'), stimulus.patch.height = 7.5;,end
if ~isfield(stimulus.patch,'sf'), stimulus.patch.sf = SF;,end
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
% function to update patch position and draw to screen
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
        
        %full screen sine contrast modulation
        sineModulation= mglMakeGrating(stimulus.patch.width,stimulus.patch.height,stimulus.patch.sf,stimulus.patch.angle,phase);
        
        %overlap annular mask
        annulus = slmakeAnnulus(stimulus.patch.width,stimulus.patch.height,stimulus.patch.inR,stimulus.patch.outR,0.5);
        
        %combine to get annulus-gabor
        gabor = convertContrast(255*(sineModulation.*annulus+1)/2);
        
        %texture for display
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
        maxX = x;
    end;
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

%draw annulus mask
function AnnulusMask = slmakeAnnulus(width,height,innerRadius,outerCircle,decayDistance)

%convert from degrees to pixels
if isempty(mglGetParam('xDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
if isempty(mglGetParam('yDeviceToPixels'))
    disp(sprintf('(makeGrating) mgl is not initialized'));
    return
end
xDeg2pix = mglGetParam('xDeviceToPixels');
yDeg2pix = mglGetParam('yDeviceToPixels');

%get annulus size in pixels
%If they are not, make width and height an odd number
widthPixels = round(width*xDeg2pix);
heightPixels = round(height*yDeg2pix);
widthPixels = widthPixels + mod(widthPixels+1,2);
heightPixels = heightPixels + mod(heightPixels+1,2);

%get a grid of x and y coordinates that has
%the correct number of pixels (grid values in degrees)
x = -width/2  : width/(widthPixels-1)   : width/2;
y = -height/2 : height/(heightPixels-1) : height/2;
[xMesh,yMesh] = meshgrid(x,y);


%%%%%%%%%%%%%%%%%%%%%%  MAKE INNER AND OUTER MASK  %%%%%%%%%%%%%%%%%%%%%%%%%
%screen center
center = [0 0];

%coord to matrix center (degrees)
x2c = xMesh - center(2);       %y to center
y2c = -(yMesh - center(1));    %x to center

%create distance-to-center map in degrees
r = sqrt(y2c.^2 + x2c.^2); %radius (Pythagoras)

%create annulus mask
AnnulusMask = ones(heightPixels,widthPixels);
AnnulusMask(r<innerRadius)=0;
AnnulusMask(r>outerCircle)=0;

%checking
checkOuter = r>outerCircle;
if ~any(checkOuter==1)
    fprintf('(slmakeAnnulus) Annulus outer edgde is too large for the screen. Reduce it.')
    mglClose; keyboard
end

%implement linear decaying to 0 over the inner .5 degree radius
%inward linear decaying mask
DecayArea = zeros(heightPixels,widthPixels);
DecayAreaIn = r>(innerRadius-decayDistance) & r<innerRadius;
DecayArea(DecayAreaIn) = r(DecayAreaIn);
x = DecayArea(DecayAreaIn);
%from out to in
scaledx = (x - min(x(:)))/(max(x(:)) - min(x(:)));
DecayArea(DecayAreaIn) = scaledx;

%outward linear decaying mask
DecayAreaOut = r<(outerCircle+decayDistance) & r>outerCircle;
DecayArea(DecayAreaOut) = r(DecayAreaOut);
x = DecayArea(DecayAreaOut);
%1 - (...) because from in to out
scaledx = 1 - ((x - min(x(:)))/(max(x(:)) - min(x(:))));
DecayArea(DecayAreaOut) = scaledx;

%Linearly decaying onward and outward Annulus
AnnulusMask = AnnulusMask + DecayArea;

%Restrict orientation angles between 1 and 180
function myangle = slKeepAngleBetween0to179(myangle)

%in deg > 0
myangle = SLra2d(SLde2r(myangle,0));
if myangle>=180
    myangle=myangle-180;
end

myangle(myangle==0) = 180;







