
% author : steeve laquitaine


function myscreen=taskDotDir(params)

%Load parameters
params = load(params);

%Initalize the screen
myscreen = initScreen;
myscreen.saveData = 1;
myscreen.datadir = cd;
myscreen.subjectID='sub01';
myscreen.screenNumber = 1;

%Hide cursor
mglDisplayCursor(0);

%Set the first task to be a fixation task
%center (in pixels)
global fixStimulus;
fixStimulus.pospix=[myscreen.screenWidth myscreen.screenHeight]/2;
[task{1},myscreen]=fixation(myscreen);

%randomize the conditions' order.
myTrialRand=randperm(numel(params.task.parameter.dir.series));

%segments of main task. Main task only have one phase.
task{2}{1}.segmin=[1 .3 5 .1 .1];
task{2}{1}.segmax=[1 .3 5 .1 .1];

%stimulus
task{2}{1}.randVars.myRandomDir=params.task.parameter.dir.series(myTrialRand);
task{2}{1}.randVars.myRandomCoh=params.task.parameter.dir.coh(myTrialRand);

%----------------------------
%for instructions (to comment)
%task{2}{1}.randVars.myRandomCoh=1;
%----------------------------

%Set initial position of the mouse randomly
task{2}{1}.randVars.initAngledeg=randi([0,359],[numel(task{2}{1}.randVars.myRandomDir),1]);

%Get response and reaction time 
%response segment is set at "3" or "4" for mouse event
task{2}{1}.getResponse=[0 0 1 0 0];
task{2}{1}.numTrials=numel(task{2}{1}.randVars.myRandomDir);

%Store data
task{2}{1}.randVars.calculated.prodcoor=[nan nan];

%Initialize our task
% task{2}{1}.thistrial.thisseg=3
[task{2}{1},myscreen]=initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);

%Init the stimulus
global stimulus;
myscreen=initStimulus('stimulus',myscreen);
stimulus=initDots(stimulus,myscreen);

%run the eye calibration
myscreen=eyeCalibDisp(myscreen);

%Main display loop
phaseNum=1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
    
    %update the dots
    [task{2},myscreen,phaseNum]=updateTask(task{2},myscreen,phaseNum);
    
    %update the fixation task
    [task{1},myscreen]=updateTask(task{1},myscreen,1);
    myscreen=tickScreen(myscreen,task);
 end

%if we got here, we are at the end of the experiment
myscreen=endTask(myscreen,task);
clear global stimulus;
clear global fixStimulus;

mglDisplayCursor(1);

%function that gets called at the start of each segment
function [task,myscreen]=startSegmentCallback(task, myscreen)
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
    %mouseinfo.buttons=mi.buttons;
    
    %if you want to you use keyboard button '1' to enter choice
    mouseinfo.buttons=mglGetKeys(19);
 
    %Position mouse on aperture
    [mouseinfo]=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);
    
    %Collect response
    task.thistrial.prodcoor=[mouseinfo.x.angle.onap mouseinfo.y.angle.onap];%(visual angle)
    [~,deg] = SLcart2polar(task.thistrial.prodcoor);
    
    %print trial number with the key and direction choosen.
    fprintf('%i %1.0f %1.0f \n',...
        task.trialnum,task.thistrial.myRandomDir,deg)
    %then jump to feedback
    task=jumpSegment(task,4);
end

%draw the stimulus each frame
function [task,myscreen]=updateScreenCallback(task, myscreen)

global stimulus;
global fixStimulus;
mglClearScreen;

%motion
if (task.thistrial.thisseg == 2)
    stimulus=updateDots(stimulus,myscreen);
end

%response
if (task.thistrial.thisseg == 3)
    
    %Position mouse on aperture
    mi=mglGetMouse(myscreen.screenNumber); %(pixels)
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    mouseinfo=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);
    
    %draw a white arrow (radius)
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        mouseinfo.x.angle.onap,...
        mouseinfo.y.angle.onap,2,[1 1 1],1);
end

%confirmation
if (task.thistrial.thisseg == 4)
    
    %Position mouse on aperture
    %deliver pixel positions, not angle
    mi=mglGetMouse(myscreen.screenNumber); 
    mouseinfo.x.pix=mi.x;
    mouseinfo.y.pix=mi.y;
    mouseinfo=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
    
    %Interface PowerMate
    [mouseinfo]=PowerMate2ap(mouseinfo,stimulus,task);   
    
    %Confirm subjects' choice by drawing a red arrow
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        mouseinfo.x.angle.onap,...
        mouseinfo.y.angle.onap,2,[1 0 0],1);
end

%feedback
if (task.thistrial.thisseg == 5)
    [coord]=polar2cartesian(task.thistrial.myRandomDir,stimulus.dots.rmax);
    feedbackcoord.x.angle.onap=coord.x;
    feedbackcoord.y.angle.onap=coord.y;
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        feedbackcoord.x.angle.onap,...
        feedbackcoord.y.angle.onap,4,[0 1 0],1);
end

%get response
function [task,myscreen]=getResponseCallBack(task,myscreen)

%When subject chose, confirm (jump to segment 4 in updateScreenCallback)
task=jumpSegment(task,4);

%init dot stimulus
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

%update dot positions and draw them to screen
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


