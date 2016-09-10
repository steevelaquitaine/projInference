%taskDotDir.m
%
%       $Id: taskDotDir.m 750 2010-03-18 02:46:46Z justin $
%     usage: taskDotDir
%        by: justin gardner modified by Steeve Laquitaine 12/12/20
%      date: 09/07/06
% copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%   purpose: random dot motion direction estimation task

%*Important basic informations*
%mglDescribeDisplays: to get information about the monitor and computer
%system
%mglDefaultScreenParams: to get info about default parameters
%mglEditScreenParams: to edit information about the monitor

%*Info about the screen in the psychophy room 121015* actually default parameters:
%- viewer distance (mm): 50.5
%- screensize(mm): 385, 285
%- reso (pix): 1280, 960
%- refresh rate: 100Hz
%*I use the following parameters for 'Marple' screen*:
%- viewer distance (mm): 50.5
%- screensize(mm): 400, 300
%- reso (pix): 960, 600: (set in setting. This resolution only allows a refresh rate of 75 Hz)

%I changed the parameters for 'Marple' screen to the 
%new parameters for the psychophysics room (121015); actually default 
%parameters.
%----------------------------------------------------------------------
%Issues
%----------------------------------------------------------------------
%[121026]
%after test task with feedbacks - it seems 90 degrees appears a lot
%more than expected... need to check it.
%cannot test the code in a small window in the main screen.
%mouse is not very comfortable to use. The laptop's pad is more
%comfortable. It takes some time to get used to the mouse.
%think about what should happen when subject doesn't respond quickly
%enough. For now I just record the last position of the mouse.
%I should only and always analyse the mouse coordinates. getangle
%function was wrong and I revised the function the 12/09/14.
%The way updatedots work. Based on my understanding of the code, the identity of the coherent dots change.
%randomly between frames. Shouldn't the same dots be updated ? Doesn't
%this procedure induce noise in the coherence ?

%[120927]
%Important issues: I register a significant number of dropped frames. 
%It artifacts the directions I perceive. The dots seem to oscillate between
%180 and 0 degrees, with a trend toward the 0).
%I tried to identify the cause of the reduced speed of the display loop.
%conclusion: 
%- mglLines2 significantly reduces the speed of the display loop and probably
%causes the drop in the number of frames that is observed.
%next steps: test the udpatescreen loop with tic toc.
%test: a whole run (internet opened) 
%dropped frames 18701 (56.78%): 42Hz refresh rate
%test: another run (internet opened)
%dropped frames 17461 (54.99%)

%[121009]
%-test: a whole run (640,480,40,30,50.5)
%dropped frames: 50%)
%-test: a whole run (960,600,40,30,50.5)
%dropped frames 12567 (53.64%): 43Hz refresh rate
%-debugging with Justin (problem solved): seg.length at 0 (randVars.myRandomCoh code) explains the dropped frames


%------------------------------------------------------------------------
%*Updates*
%------------------------------------------------------------------------
%[121015] later task updates:
%cleaner data with coh: 100%, 50%, 24%, 8%.

%[121026]
%%feedbacks about true direction added
%%arrows are antialiased because it can be used to reproduced the
%direction as a benchmark.i.e., lines are perfectly linear for the
%cardinals and obliques multiples of 45 degrees.

%[121201]
%Initial position of response arrow is set at (prior mean + cardinal)/2. 
%In this way the effect of initial position is counterbalanced between prior
%cardinal. A disadvantage is that it can weaken both the attraction of the prior and
%the cardinal. But importantly, its effect is cancelled when contrasting prior
%conditions. When using 2AFC, we won't have to worry about such bias.

%[121202]
%- I should use a coherence between 0.24 and 1 to observe a good qualitative
%effect of the conflict between prior (upward from the diagonal) and
%cardinal(downward from the diagonal);
%!!!---I get the max expected effect at coh: 0.35---!!!.
%- I think coh: 0.08 is too low. We should use at minimum 0.24. I perceive
%strange motion effects (e.g., motion direction waves).

%[130209]
%faster trial: 
%fixation: 0.7 s
%motion: 0.35 s
%feedback: 0.5 s
%powermate motion faster: wheelspeed 0.016 instead of 0.02.
%feeback arrow width doubled: 4 instead of 2: to reinforce the prior

%[130806]
%get reaction time
%when subject doesn't respond within time limits, reaction time marks NaN.

%[131220]
%running this code without any change produced error message probably due
%to an update in an mgl function. I correct the code.
%updated taskDotDir code now set directory to save file and subjectID


function myscreen=taskDotDir(params,saveinfo)

%Check arguments
if ~any(nargin==1)
    help taskDotDir
    return
end

%Check if parameters exist in directory
if isempty(dir([params,'.mat']))==1; 
disp([{' ---- your parameters were not found in the directory ---- '}; 
    {'does "taskDotDir" varargin corresponds to a file in the directory ?'}])    
return; 
end

%Load parameters
params=load(params);

%Initalize the screen
%for debugging: uncomment
%mglOpen(0)
%mglMoveWindow(100,800)
myscreen=initScreen;
myscreen.saveData=1;
myscreen.datadir=saveinfo.datadir;
myscreen.subjectID=saveinfo.subjectID;

%Set the target screen
myscreen.screenNumber=2;

%Hide cursor
mglDisplayCursor(0);

%Set the first task to be a fixation task
global fixStimulus;
fixStimulus.pospix=[myscreen.screenWidth myscreen.screenHeight]/2; %center position (pixels)
[task{1},myscreen]=fixation(myscreen);

%Set the second task to have five segments.
%task{2}{1}.segmin=[.5 .35 5 .15 .2]; %motion and feedback 700 ms
%task{2}{1}.segmax=[.5 .35 5 .15 .2];
task{2}{1}.segmin=[1 .3 5 .1 .1]; %motion and feedback 700 ms
task{2}{1}.segmax=[1 .3 5 .1 .1];

%task{2}{1}.parameter.dir=params.task.parameter.dir.series;%repmat(130,1,65) %time series comes from a .mat
%task{2}{1}.parameter.coherence=[.12 .35 1]; %new pair of 3 coherences are tested (clearer effects).
%Randomize the conditions' order.
myTrialRand=randperm(numel(params.task.parameter.dir.series));
%Extract the stimulus properties.
task{2}{1}.randVars.myRandomDir=params.task.parameter.dir.series(myTrialRand);%repmat(130,1,65) %time series comes from a .mat
task{2}{1}.randVars.myRandomCoh=params.task.parameter.dir.coh(myTrialRand);%[.12 .35 1]; 3 coherences are tested (clearer effects).
%task{2}{1}.random=1;

%Get response and reaction time 
%response segment is set at "4".
task{2}{1}.getResponse=[0 0 4 0 0];

%task{2}{1}.numTrials=numel(task{2}{1}.randVars.myRandomDir)*...
%numel(task{2}{1}.randVars.myRandomCoh); %all trials is ideal choice.
task{2}{1}.numTrials=numel(task{2}{1}.randVars.myRandomDir);%all trials

%Set initial position of the mouse randomly
task{2}{1}.randVars.initAngledeg=randi([0,359],[numel(task{2}{1}.randVars.myRandomDir),1]);

%Store behavioral data
task{2}{1}.randVars.calculated.prodcoor=[nan nan];

%Initialize our task
for phaseNum=1:length(task{2})
    [task{2}{phaseNum},myscreen]=initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);
end

%Init the stimulus
global stimulus;
myscreen=initStimulus('stimulus',myscreen);
stimulus=initDots(stimulus,myscreen);
%%debugging code - remove later
%%allocate space for 10 minutes worth of data for 5 different checkpoints
%stimulus.tictoc= nan(8,60*60*10);
%stimulus.tictoci=1;
%stimulus.tictoci2=1;
%stimulus.tictoci3=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen=eyeCalibDisp(myscreen);

%-------------------------------------------------------------------
%Main display loop
%-------------------------------------------------------------------
phaseNum=1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
    %update the dots
    %startTime=mglGetSecs;
    [task{2},myscreen,phaseNum]=updateTask(task{2},myscreen,phaseNum);
%     [task{2},myscreen,phaseNum]=updateTask(task{2},myscreen,phaseNum);
    %stimulus.tictoc(6,stimulus.tictoci3)=mglGetSecs(startTime);
    %update the fixation task
    [task{1},myscreen]=updateTask(task{1},myscreen,1);
    %stimulus.tictoc(7,stimulus.tictoci3)=mglGetSecs(startTime);
    %flip screen
    myscreen=tickScreen(myscreen,task);
    %stimulus.tictoc(8,stimulus.tictoci3)=mglGetSecs(startTime);
    %stimulus.tictoci3=stimulus.tictoci3+1;
 end

%if we got here, we are at the end of the experiment
myscreen=endTask(myscreen,task);

%%debugging 
%mglDisplayCursor(1);
%for i=1:size(stimulus.tictoc,1)
%  disp(sprintf('Mean time for %i: %0.4f ms Max: %0.4f',i,1000*nanmean(stimulus.tictoc(i,:)),1000*nanmax(stimulus.tictoc(i,:))));
%end
%
%keyboard

clear global stimulus;
clear global fixStimulus;

%function that gets called at the start of each segment
function [task,myscreen]=startSegmentCallback(task, myscreen)
%startTime=mglGetSecs;

global stimulus;
global fixStimulus;
global Inicoord;
%motion
if (task.thistrial.thisseg == 2)
    stimulus.dots.coherence=task.thistrial.myRandomCoh;
end
stimulus.dots.dir=task.thistrial.myRandomDir;

%response
if (task.thistrial.thisseg == 3)
    %Position mouse in the center of the screen
    %mglSetMousePosition(fixStimulus.pospix(1), fixStimulus.pospix(2), myscreen.screenNumber);
    %Calculate angle coordinates on aperture relative to center
    [coord]=polar2cartesian((225+270)/2,stimulus.dots.rmax);
    Inicoord.x.angle.onap=coord.x;
    Inicoord.y.angle.onap=coord.y;
    %Calculate pixels coordinates on aperture relative to center
    Inicoord.x.pix.onap=Inicoord.x.angle.onap * mglGetParam('xDeviceToPixels');
    Inicoord.y.pix.onap=Inicoord.y.angle.onap * mglGetParam('yDeviceToPixels');
    %Calculate coordinates on aperture relative to screen's root.
    Inicoord2root.x.pix.onap=Inicoord.x.pix.onap + fixStimulus.pospix(1);
    Inicoord2root.y.pix.onap=Inicoord.y.pix.onap + fixStimulus.pospix(2);
    %Position mouse
    mglSetMousePosition(Inicoord2root.x.pix.onap, Inicoord2root.y.pix.onap, myscreen.screenNumber);   
end

%confirmation    
if (task.thistrial.thisseg == 4)
    %Get mouse position
    mi=mglGetMouse(myscreen.screenNumber); %get pixel positions.  
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    
    %Check if subject confirmed his choice: ("space bar is down")
    %if you want to you use mouse click to enter choice
    mouseinfo.buttons=mi.buttons;%mglGetKeys(50); %mi.buttons;
%    %if you want to you use keyborad button '0' to enter choice
%    mouseinfo.buttons=mglGetKeys(83);%mglGetKeys(50) %mi.buttons; 
%    
    %Position mouse on aperture
    [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);
    %Collect response
    task.thistrial.prodcoor=[mouseinfo.x.angle.onap mouseinfo.y.angle.onap];%(visual angle)
end
%It would be nice if i could implement a static feedback here instead of
%UpdateScreenCallback function, but I don't understand how to use mglflush.mode
%function.

%function that gets called to draw the stimulus each frame
function [task,myscreen]=updateScreenCallback(task, myscreen)

global stimulus;
global fixStimulus;

%startTime=mglGetSecs;
%global mouseinfo;
mglClearScreen;
%stimulus.tictoc(1,stimulus.tictoci)=mglGetSecs(startTime);

%motion
if (task.thistrial.thisseg == 2)
    stimulus=updateDots(stimulus,myscreen);
end
%stimulus.tictoc(2,stimulus.tictoci)=mglGetSecs(startTime);

%response
if (task.thistrial.thisseg == 3)
    %translate mouse position 2 the aperture
    mi=mglGetMouse(myscreen.screenNumber); %(pixels)
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;    
    %Check if subject confirmed his choice: ("space bar is down")
    %if you want to you use mouse click to enter choice
    mouseinfo.buttons=mi.buttons;%mglGetKeys(50); %mi.buttons;
%  %if you want to you use keyborad button '0' to enter choice
%    mouseinfo.buttons=mglGetKeys(83); %mi.buttons;
    
    %Position mouse on aperture
    [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);
    %when no response
    if ~mouseinfo.buttons     
        %draw a white arrow linking the center to the mouse position
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            mouseinfo.x.angle.onap,...
            mouseinfo.y.angle.onap,2,[1 1 1],1);
        %at response, confirm response.
    else
        task=jumpSegment(task,4);
    end
end
%0.25 ms/frames at max here
%stimulus.tictoc(3,stimulus.tictoci)=mglGetSecs(startTime);

%confirmation
if (task.thistrial.thisseg == 4) %confirmation (segment 4)
    %translate mouse 2 aperture
    mi=mglGetMouse(myscreen.screenNumber); %deliver pixel positions, not angle    
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    
    %if you want to you use mouse click to enter choice
    mouseinfo.buttons=mi.buttons;%mglGetKeys(50); %mi.buttons;
%    %if you want to you use keyborad button '0' to enter choice
%    mouseinfo.buttons=mglGetKeys(83); %mi.buttons;
    
    %Position mouse on aperture
    [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);   
    %Confirm subjects' choice by drawing a red arrow
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        mouseinfo.x.angle.onap,...
        mouseinfo.y.angle.onap,2,[1 0 0],1);
end
%debugging
%stimulus.tictoc(4,stimulus.tictoci)=mglGetSecs(startTime);
%stimulus.tictoci=stimulus.tictoci+1;
%0.29 ms/frames at max here
%0.79 ms/frames at max here

%feedback
if (task.thistrial.thisseg == 5)
    %debugging - to remove
%    task.thistrial.myRandomDir=90;  
    [coord]=polar2cartesian(task.thistrial.myRandomDir,stimulus.dots.rmax);
    feedbackcoord.x.angle.onap=coord.x;
    feedbackcoord.y.angle.onap=coord.y;
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        feedbackcoord.x.angle.onap,...
        feedbackcoord.y.angle.onap,4,[0 1 0],1);
%    disp(task.thistrial.myRandomDir)
end

%function that gets called to get response
function [task,myscreen]=getResponseCallBack(task,myscreen)
task.thistrial.buttonState

%function to init the dot stimulus
function stimulus=initDots(stimulus,myscreen)

%convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax=2.5;end %(radius in angle)
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter=0;end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter=0;end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize=3;end %(in pixels)
if ~isfield(stimulus.dots,'density'), stimulus.dots.density=16.7;end %dots/degree/sec (Hanks et al,2012,JN)
%if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed=1.5;,end %2 (previous exp.)%8;,end %Rao et al, 2012,JN
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed=2.8;end %2 (previous exp.)%8;,end %Rao et al, 2012,JN


%actually a square patch of dots that get stenciled
%so calculate width and height
stimulus.dots.width=stimulus.dots.rmax*2;
stimulus.dots.height=stimulus.dots.rmax*2;

%get the number of dots
stimulus.dots.n=round(stimulus.dots.width*stimulus.dots.height*stimulus.dots.density);

%get max and min points for dots
stimulus.dots.xmin=-stimulus.dots.width/2;
stimulus.dots.xmax=stimulus.dots.width/2;
stimulus.dots.ymin=-stimulus.dots.height/2;
stimulus.dots.ymax=stimulus.dots.height/2;

%get initial position
stimulus.dots.x=rand(1,stimulus.dots.n)*stimulus.dots.width;
stimulus.dots.y=rand(1,stimulus.dots.n)*stimulus.dots.height;

%get the step size
stimulus.dots.stepsize=stimulus.dots.speed/myscreen.framesPerSecond;

%create stencil
mglClearScreen;
mglStencilCreateBegin(1);
%and draw that oval
%mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax]/2,[1 1 1],60);
mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax],[1 1 1],60);
mglStencilCreateEnd;
mglClearScreen;

%function to update dot positions and draw them to screen
function stimulus=updateDots(stimulus,myscreen)

%get the dots step
stimulus.dots.xstep=cos(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;
stimulus.dots.ystep=sin(pi*stimulus.dots.dir/180)*stimulus.dots.stepsize;

%pick a random set of dots
stimulus.dots.coherent=rand(1,stimulus.dots.n) < stimulus.dots.coherence;

%now move those dots in the right direction
stimulus.dots.x(stimulus.dots.coherent)=stimulus.dots.x(stimulus.dots.coherent)+stimulus.dots.xstep;
stimulus.dots.y(stimulus.dots.coherent)=stimulus.dots.y(stimulus.dots.coherent)+stimulus.dots.ystep;

%randomwalk rule
thisdir=rand(1,sum(~stimulus.dots.coherent))*2*pi;
stimulus.dots.x(~stimulus.dots.coherent)=stimulus.dots.x(~stimulus.dots.coherent)+cos(thisdir)*stimulus.dots.stepsize;
stimulus.dots.y(~stimulus.dots.coherent)=stimulus.dots.y(~stimulus.dots.coherent)+sin(thisdir)*stimulus.dots.stepsize;

%movshon noise
%stimulus.dots.x(~stimulus.dots.coherent)=rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.width;
%stimulus.dots.y(~stimulus.dots.coherent)=rand(1,sum(~stimulus.dots.coherent))*stimulus.dots.height;

%make sure we haven't gone off the patch
%do the dots separately for left and right hand side
stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)=stimulus.dots.x(stimulus.dots.x < stimulus.dots.xmin)+stimulus.dots.width;
stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)=stimulus.dots.x(stimulus.dots.x > stimulus.dots.xmax)-stimulus.dots.width;
stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)=stimulus.dots.y(stimulus.dots.y < stimulus.dots.ymin)+stimulus.dots.height;
stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)=stimulus.dots.y(stimulus.dots.y > stimulus.dots.ymax)-stimulus.dots.height;

%draw the dots
mglStencilSelect(1);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,[1 1 1]);
mglStencilSelect(0);





%function to convert from mouse's to aperture's coordinates
function [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus)

%and the mouse coordinates relative to the center of the target screen
mouseinfo.x.pix=mouseinfo.x.pix-fixStimulus.pospix(1); %(in pix)
mouseinfo.y.pix=mouseinfo.y.pix-fixStimulus.pospix(2);

%convert pixels 2 angles (mglLines works with angles)
mouseinfo.x.angle.screen=mouseinfo.x.pix*mglGetParam('xPixelsToDevice');
mouseinfo.y.angle.screen=mouseinfo.y.pix*mglGetParam('yPixelsToDevice');

%calculate the transformation parameter that position the mouse on the aperture
%calculate length of the arrow (pythagoras theorem)
arrowsz2=mouseinfo.x.angle.screen^2 + mouseinfo.y.angle.screen^2;%(angle)

%calculate transformation parameter
transpara=stimulus.dots.rmax^2/arrowsz2;

%transform actual coordinates to on-aperture coordinates.
mouseinfo.x.angle.onap=sqrt(mouseinfo.x.angle.screen^2*transpara)*sign(mouseinfo.x.angle.screen);
mouseinfo.y.angle.onap=sqrt(mouseinfo.y.angle.screen^2*transpara)*sign(mouseinfo.y.angle.screen);
%function to convert from PowerMate's to aperture's coordinates
function [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task)
global Inicoord;
%Set speed, initial position and size of arrow.
wheelspeed=0.008;
%fixed
%initAngle=(225+270)/2*pi/180;
%random
initAngle.rad=task.thistrial.initAngledeg*pi/180;

r=stimulus.dots.rmax;
%Calculate arrow's coordinates on aperture
mouseinfo.x.angle.onap= r * cos(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
mouseinfo.y.angle.onap= r * sin(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap));
%function to convert from polar to cartesian coordinates
function [coord]=polar2cartesian(theta,r)
%theta is an angle in degree
%r is the radius of the unit circle
%Coord are in visual angle
%Record angle in degree
theta2.deg=theta;
%Convert from degree to radian
theta2.rad=theta2.deg*pi/180;
%Calculate visual angles coordinates
coord.x=r*cos(theta2.rad);
coord.y=r*sin(theta2.rad);
