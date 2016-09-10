%SLfmriLocalizerMotion.m
%
%       $Id: SLfmriLocalizerMotion.m,v 1.6 2009/02/02 03:55:00 justin Exp $
%        by: steeve laquitaine based on Justin Gardner's codes, (c) 2006 
%            Justin Gardner (GPL see mgl/COPYING)
%      date: 09/07/06 updated 150112
%   purpose: Motion localizer fMRI experiment
%
%     usage:
%
%           SLfmriLocalizerMotion('randMotion')
%           SLfmriLocalizerMotion('black')
%
%Description
%
%   waitForBacktick option is set to 1 which 
%
%   task{2}{1}.waitForBacktick = 1;

function myscreen = SLfmriLocalizerMotion(type)

%data folder
myscreen.datadir = cd;

%initalize black screen
myscreen.autoCloseScreen = 1;
myscreen.saveData=1; 
myscreen.displayname='projector';
myscreen.background='black';
myscreen.screenNumber=1;
myscreen=initScreen(myscreen);

%set the first task to be the fixation staircase task
%When it's a staircase the fixation is surrounded by a black circular patch
%of 3? diameter (KoK et al., 2013)
%The fixation task is a 2AFC fixation dimming task in which the fixation 
%cross dims twice and the subject has to respond with either a button press
%on 1 or 2 to indicate which interval was dimmer. Subjects answer when the
%fixation turn yellow and is instructed of correct choice (yellow turn 
%green) or incorrect choice (yellow turn red).
%This is run on a staircase so that the difficulty is adjusted throughout 
%a scan to maintain vigilance and attention level.
%run this code for instructed training: "fixStairTest".
clear global fixStimulus;
global fixStimulus;
fixStimulus.diskSize=1.5;
%[task{1},myscreen]=fixation(myscreen);
[task{1} myscreen] = fixStairInitTask(myscreen);

%a top-up period of the same direction
%The localizer is composed of 11 trials of 25.12s (4min 36s); 176 backticks
%every 1.57s(TR);
%-each cycle(trial) is composed of 16 volumes of 1.57s (TR).
%-each cycle(trial) is composed of 6 segments (5 coherent motion,1 random)
%-the first trial (16 volumes) is discarded (period of BOLD stabilisation).
%Volume acquisition is re-synchronized with the task after the segment 
%for which synchToVol=1. The backtick signal is sent at each volume 
%(every 1.57s (TR)) and is resynchronized at the start of every cycle.
%So synchToVol must be set at 1 before the end of trial's last volume. 
%The last volume starts ends at 25.12s, so a good trial duration is
%24s.
%to simulate the task during fMRI acquisition run 
%mglSimulateRun(1.57,176,10); 
%dotslocnewSteeve

%fMRI acquisition should finish exactly at 11 trials but we add 10 more
%trials (50.24s) to be sure we get everything.
%How should I choose trial duration given a number of vols/trial, and that
%TR ranges between 1.5 & >1.5 ?
%The number of volumes must be chosen correctly relative to trial duration, 
%e.g., If I want to acquire 11 trials with 176 volumes (16 vols/trials) given 
%that the TR changes such that 1.5s<=TR<=1.71s, a trial should last 
%24s (because 24/16 vols <= TR <= 24/15 vols).
%e.g., 176 vols with a short 1.5s TR doesn't permit to aquire 11 trials of
%25s.
task{2}{1}.seglen     = [2.4 2.4 2.4 2.4 2.4 12];
task{2}{1}.synchToVol = [  0   0   0   0   0  1];
task{2}{1}.parameter.coherence = 1;
task{2}{1}.random = 1;
task{2}{1}.numTrials = 11 + 10;
task{2}{1}.fudgeLastVolume = 1;
task{2}{1}.waitForBacktick = 1;

%initialize our task
for phaseNum=1:length(task{2})
  [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback);
end

%init the stimulus
clear global stimulus;
global stimulus;
myscreen=initStimulus('stimulus',myscreen);
stimulus.dots.type='Opticflow';
stimulus=initDots(stimulus,myscreen);
stimulus.type=type;

%run the eye calibration
myscreen = eyeCalibDisp(myscreen);


%Main display loop
phaseNum=1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
  %update the dots
  [task{2} myscreen phaseNum]=updateTask(task{2},myscreen,phaseNum);
  %update the fixation task
  [task{1} myscreen]=updateTask(task{1},myscreen,1);
  %flip screen
  myscreen=tickScreen(myscreen,task);
end

%if we got here, we are at the end of the experiment
myscreen=endTask(myscreen,task);


%function that gets called at the start of each segment
function [task myscreen]=startSegmentCallback(task, myscreen)
global stimulus;
stimulus.dots.color=[1 1 1];

%when the segment is not the last segment
if (task.thistrial.thisseg < length(task.seglen))
  stimulus.coherence=task.thistrial.coherence;
  %set speed
  stimulus.dots=feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,stimulus.speed,myscreen);
  stimulus.dots=feval(sprintf('setDotsDir%s',stimulus.dots.type),stimulus.dots,2*mod(task.thistrial.thisseg,2)-1,myscreen);
else    
    if strcmp(stimulus.type,'randMotion')
        stimulus.coherence=0;
    end
    if strcmp(stimulus.type,'static')
        stimulus.coherence=0;
        stimulus.dots=feval(sprintf('setDotsSpeed%s',stimulus.dots.type),stimulus.dots,0,myscreen);
    end
    if strcmp(stimulus.type,'black')
        stimulus.coherence=0;
        stimulus.dots.color=[0 0 0];
    end
end

%function that gets called to draw the stimulus each frame
function [task myscreen]=updateScreenCallback(task, myscreen)

global stimulus
mglClearScreen;

%update the dots
stimulus.dots=feval(sprintf('updateDots%s',stimulus.dots.type),stimulus.dots,stimulus.coherence,myscreen);

%draw the dots
if stimulus.dots.mask,mglStencilSelect(1);end
%mglPoints2(stimulus.dots.x(stimulus.dots.color==1),stimulus.dots.y(stimulus.dots.color==1),stimulus.dots.dotsize,[1 1 1]);
%mglPoints2(stimulus.dots.x(stimulus.dots.color==0),stimulus.dots.y(stimulus.dots.color==0),stimulus.dots.dotsize,[0 0 0]);
mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,stimulus.dots.color);
if stimulus.dots.mask,mglStencilSelect(0);end

%mask the dots for the last segment (black screen)
if task.thistrial.thisseg == length(task.seglen)
    if stimulus.dots.mask,mglStencilSelect(1);end
    mglPoints2(stimulus.dots.x,stimulus.dots.y,stimulus.dots.dotsize,stimulus.dots.color);
    if stimulus.dots.mask,mglStencilSelect(0);end
end

%function to init the dot stimulus
function stimulus=initDots(stimulus,myscreen)

%convert the passed in parameters to real units
% if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax=2.5;
% if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter=0;end
% if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter=0;end
% if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize=3;end %(in pixels)
% if ~isfield(stimulus.dots,'density'), stimulus.dots.density=16.7;end %dots/degree^2/frames (Hanks et al,2012,JN)
% %if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed=1.5;,end %2 (previous exp.)%8;,end %Rao et al, 2012,JN
% if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed=2.8;end %2 (previous exp.)%8;,end %Rao et al, 2012,JN
%fMRI
%based on a mixture of Kok et al,JN,2013 and psychophysics task.
%the diameter is 15? as in Kok et al,JN,2013 and Vintch et al,
%JN,2013 (0.5-20?radius? probably diameter)
%psychophysics
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax=11;end;
if ~isfield(stimulus.dots,'xcenter'), stimulus.dots.xcenter=0;,end
if ~isfield(stimulus.dots,'ycenter'), stimulus.dots.ycenter=0;,end
if ~isfield(stimulus.dots,'dotsize'), stimulus.dots.dotsize=3;,end
if ~isfield(stimulus.dots,'density'), stimulus.dots.density=16.7;,end
if ~isfield(stimulus.dots,'coherence'), stimulus.dots.coherence=1;,end
stimulus.speed=6;
if ~isfield(stimulus.dots,'speed'), stimulus.dots.speed=stimulus.speed;,end
if ~isfield(stimulus.dots,'dir'), stimulus.dots.dir=0;,end
if ~isfield(stimulus.dots,'mask'), stimulus.dots.mask=1;,end

%update the dots
stimulus.dots=feval(sprintf('initDots%s',stimulus.dots.type),stimulus.dots,myscreen);

%set color
stimulus.dots.color=ones(stimulus.dots.n,1);
%stimulus.dots.color(rand(1,stimulus.dots.n)>0.5)=1;

%create stencil
if stimulus.dots.mask
  mglClearScreen;
  mglStencilCreateBegin(1);
  %and draw that oval
  mglGluDisk(stimulus.dots.xcenter,stimulus.dots.ycenter,[stimulus.dots.rmax stimulus.dots.rmax],[1 1 1],60);
  mglStencilCreateEnd;
  mglClearScreen;
end

%set the dots speed
function dots=setDotsSpeedLinear(dots,speed,myscreen)

%get the step size
dots.speed=speed;
dots.stepsize=speed/myscreen.framesPerSecond;

%set the dots direction
function dots=setDotsDirectionLinear(dots,direction,myscreen)

%get the step size
dots.dir=direction;

%step dots for linear
function dots=updateDotsLinear(dots,coherence,myscreen)

%get the dots step
dots.xstep=cos(pi*dots.dir/180)*dots.stepsize;
dots.ystep=sin(pi*dots.dir/180)*dots.stepsize;

%pick a random set of dots
dots.coherent=rand(1,dots.n) < coherence;

%now move those dots in the right direction
dots.x(dots.coherent)=dots.x(dots.coherent)+dots.xstep;
dots.y(dots.coherent)=dots.y(dots.coherent)+dots.ystep;

%randomwalk rule
thisdir=rand(1,sum(~dots.coherent))*2*pi;
dots.x(~dots.coherent)=dots.x(~dots.coherent)+cos(thisdir)*dots.stepsize;
dots.y(~dots.coherent)=dots.y(~dots.coherent)+sin(thisdir)*dots.stepsize;

%movshon noise
%dots.x(~dots.coherent)=rand(1,sum(~dots.coherent))*dots.width;
%dots.y(~dots.coherent)=rand(1,sum(~dots.coherent))*dots.height;

%make sure we haven't gone off the patch
%do the dots separately for left and right hand side
dots.x(dots.x < dots.xmin)=dots.x(dots.x < dots.xmin)+dots.width;
dots.x(dots.x > dots.xmax)=dots.x(dots.x > dots.xmax)-dots.width;
dots.y(dots.y < dots.ymin)=dots.y(dots.y < dots.ymin)+dots.height;
dots.y(dots.y > dots.ymax)=dots.y(dots.y > dots.ymax)-dots.height;

%step dots for opticflow
function dots=updateDotsOpticflow(dots,coherence,myscreen)

%get the coherent and incoherent dots
%if (dots.coherency ~= coherence)
  dots.incoherent=rand(1,dots.n) > coherence;
  dots.incoherentn=sum(dots.incoherent);
  dots.coherent=~dots.incoherent;
  dots.coherency=coherence;
  %generate a random transformation matrix for each incoherent point
  dots.randT=rand(3,dots.incoherentn)-0.5;
  %and normalize the transformation to have the same length
  %(i.e. speed) as the real transformation matrix
  dots.randT=sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
%end

%update relative position of dots in 3-space to observer
dots.X(dots.coherent)=dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent)=dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent)=dots.Z(dots.coherent)-dots.T(3);

%now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent)=dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent)=dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent)=dots.Z(dots.incoherent)-dots.randT(3,:);

%get all points that have fallen off the screen
offscreen=dots.Z<dots.minZ;

%and put them at the furthest distance
dots.Z(offscreen)=dots.maxZ;

%get all points that have fallen out of view
offscreen=dots.Z>dots.maxZ;
%and move them to the front plane
dots.Z(offscreen)=dots.minZ;

%put points fallen off the X edge back
offscreen=dots.X < -dots.maxX;
dots.X(offscreen)=dots.X(offscreen)+2*dots.maxX;
offscreen=dots.X > dots.maxX;
dots.X(offscreen)=dots.X(offscreen)-2*dots.maxX;

%put points fallen off the Y edge back
offscreen=dots.Y < -dots.maxY;
dots.Y(offscreen)=dots.Y(offscreen)+2*dots.maxY;
offscreen=dots.Y > dots.maxY;
dots.Y(offscreen)=dots.Y(offscreen)-2*dots.maxY;

%project on to screen
dots.xproj=dots.f*dots.X./dots.Z;
dots.yproj=dots.f*dots.Y./dots.Z;

%stuff to compute median speed
dots.oldx=dots.x;
dots.oldy=dots.y;

%get actual screen coordinates
dots.x=dots.xproj*myscreen.imageWidth;
dots.y=dots.yproj*myscreen.imageHeight;

%medianSpeed=median(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%minSpeed=min(sqrt((dots.oldx-dots.x).^2+(dots.oldy-dots.y).^2)*myscreen.framesPerSecond);
%disp(sprintf('min: %f median: %f',minSpeed,medianSpeed));

%create dots for linear2
function dots=initDotsLinear(dots,myscreen)

%actually a square patch of dots that get stenciled
%so calculate width and height
dots.width=dots.rmax*2;
dots.height=dots.rmax*2;

%get the number of dots
dots.n=round(dots.width*dots.height*dots.density);

%get max and min points for dots
dots.xmin=-dots.width/2;
dots.xmax=dots.width/2;
dots.ymin=-dots.height/2;
dots.ymax=dots.height/2;

%get initial position
dots.x=rand(1,dots.n)*dots.width;
dots.y=rand(1,dots.n)*dots.height;

%get the step size
dots.stepsize=dots.speed/myscreen.framesPerSecond;

%set the dots speed
function dots=setDotsSpeedOpticflow(dots,speed,myscreen)

%get the step size
dots.speed=speed;
dots.T=[0 0 dots.speed/myscreen.framesPerSecond];

%set the dots direction
function dots=setDotsDirOpticflow(dots,direction,myscreen)

%get the step size
dots.T=[0 0 direction*dots.speed/myscreen.framesPerSecond];

%create dots for optic flow
function dots=initDotsOpticflow(dots,myscreen)

%focal length to projection plane
%projection plane is defined to be 
%1 unit wide and high, so with 
%this focal length, we are looking at
%a view of the world with a 90 deg fov
dots.f=.5;

%translation and rotation matrices
dots.T=[0 0 dots.speed/myscreen.framesPerSecond];
dots.R=[0 0 0];

%maximum depth of points
dots.maxZ=10;dots.minZ=dots.f;
dots.maxX=10;
dots.maxY=10;

%make a brick of points
dots.n=round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

%initial position of dots
dots.X=2*dots.maxX*rand(1,dots.n)-dots.maxX;
dots.Y=2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z=(dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

%get projection on to plane
dots.xproj=dots.f*dots.X./dots.Z;
dots.yproj=dots.f*dots.Y./dots.Z;

%put into screen coordinates
dots.x=dots.xproj*myscreen.imageWidth;
dots.y=dots.yproj*myscreen.imageHeight;

%set incoherent dots to 0
dots.coherency=1;
dots.incoherent=rand(1,dots.n) > dots.coherency;
dots.incoherentn=sum(dots.incoherent);
dots.coherent=~dots.incoherent;

dots.randT=zeros(3,dots.incoherentn);