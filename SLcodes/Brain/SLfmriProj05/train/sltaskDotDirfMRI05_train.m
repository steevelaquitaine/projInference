%sltaskDotDirfMRI05_train.m
%note: oban
%        by: Steeve Laquitaine inspired by justin gardner codes
%      date: 2016/04/23
%   purpose: random dot motion direction estimation task to train
%            subject with one prior distribution of motion directions 
%            before each fMRI session
%     usage:
%
%           sltaskDotDirfMRI05_train('fMRIprojFlex','TRAINparamsPriors80m225Coh12and06',1,3)
%
%   inputs:
%       RespConfirm : 1 (Mouse click or CNI touch pad click)
%                   : 2 (keyboard button '1')
%
%          RespType : how to move the response arrow
%                     - 3 is for keyboard buttons 1 (ccw) and 2 (cw)
%
%Instructions:
%     - You will view a field of dots moving more or less coherently in a given direction.
%     - When a response arrow appears indicate the coherent direction of the motion by moving
%       the mouse horizontally to adjust the arrow position 
%       (adjust instructions according to choosen response type).
%     - Click the mouse button (if RespConfirm=1) or press keyboard key 1 (if RespConfirm=2)
%       or press "3" (if RespConfirm=3) to enter your estimate.
%     - a red arrow will confirm the direction of your choice
%     - followed by a green arrow that indicates the true direction of the
%       motion
%     - Repeat for each trial        


function myscreen = sltaskDotDirfMRI05_train(screenName,paramsFile,RespConfirm,RespType)

%Load parameters
params = load(paramsFile);
fprintf('%s %s \n','(slTaskDotDirfMRI05) Loading ',paramsFile)

%Initialize screen
myscreen = initScreen(screenName);
myscreen.screenNumber(myscreen.screenNumber==0)=1;
fprintf('%s \n','(sltaskDotDirfMRI05) screenNumber 0 re-set to 1 because testing')
myscreen.saveData = 1;
myscreen.RespConfirm = RespConfirm;  %1 (Mouse or CNI touchpad) 2 (Powermate)
myscreen.responseType = RespType;
myscreen.RespConfirm(myscreen.responseType==3)=3;
mglDisplayCursor(1); %Hide cursor
myscreen.keyboard.nums = [19 20 21]; %listen to keys 1,2,3

%Set the first task to be a fixation task
global fixStimulus;
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2; %center position (pixels)
[task{1},myscreen] = fixationTrain(myscreen);

%Set the second task to have five segments.
task{2}{1}.segmin = [1 .3 5 .1 .1]; %motion and feedback 700 ms
task{2}{1}.segmax = [1 .3 5 .1 .1];

%Randomize the conditions' order.
myTrialRand = randperm(numel(params.task.parameter.dir.series));
%dot task params
task{2}{1}.randVars.myRandomDir  = params.task.parameter.dir.series(myTrialRand);

%----------------
%adjust coherence
%----------------
params.task.parameter.dir.coh(params.task.parameter.dir.coh==0.12) = 0.24;
params.task.parameter.dir.coh(params.task.parameter.dir.coh==0.06) = 0.12;
task{2}{1}.randVars.myRandomCoh  = params.task.parameter.dir.coh(myTrialRand);


%task{2}{1}.randVars.myRandomCoh  = ones(1,length(myTrialRand)); %forpractice
task{2}{1}.randVars.myMean       = params.task.parameter.dir.mean;    %prior mean
task{2}{1}.randVars.myModes      = params.task.parameter.dir.mean;    %prior mode
task{2}{1}.randVars.myStrength   = 80;                                %prior strength
task{2}{1}.randVars.myDirections = params.task.parameter.dir.sample;  %unique directions

%Get response and reaction time
%response segment is set at "3" or 4" for mouse event
%response segment is set at "1" for keyboard event
if myscreen.RespConfirm == 1
    task{2}{1}.getResponse = [0 0 4 0 0];
    fprintf('(slTaskDotDirfMRI05_train) Response confirmation is set with Mouse or CNI touchpad click  \n')
elseif myscreen.RespConfirm == 2
    task{2}{1}.getResponse = [0 0 1 0 0];
    fprintf('(slTaskDotDirfMRI05_train) Response confirmation is set with keyboard key "1" \n')
elseif myscreen.RespConfirm == 3
    task{2}{1}.getResponse = [0 0 1 0 0];
    fprintf('(slTaskDotDirfMRI05_train) Response confirmation is set with keyboard key "3" \n')
end

%task{2}{1}.numTrials=numel(task{2}{1}.randVars.myRandomDir)*...
%numel(task{2}{1}.randVars.myRandomCoh); %all trials is ideal choice.
task{2}{1}.numTrials = numel(task{2}{1}.randVars.myRandomDir);%all trials

%Set initial position of the mouse randomly
task{2}{1}.randVars.initAngledeg = randi([0,359],[numel(task{2}{1}.randVars.myRandomDir),1]);

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

%response arrow length (deg)
stimulus.rmaxResp = 1;
stimulus.respSpeed = 0.1;

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
function [task,myscreen] = startSegmentCallback(task, myscreen)

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
    
    %initialize confirmation button state
    task.thistrial.buttonConfstate = 0;
    
    %Position mouse in the center of the screen
    mglSetMousePosition(fixStimulus.pospix(1), fixStimulus.pospix(2), myscreen.screenNumber);
    
    %Calculate angle coordinates on aperture relative to center
    [coord] = polar2cartesian((225+270)/2,stimulus.dots.rmax);
    Inicoord.x.angle.onap = coord.x;
    Inicoord.y.angle.onap = coord.y;
    
    %Calculate pixels coordinates on aperture relative to center
    Inicoord.x.pix.onap=Inicoord.x.angle.onap * mglGetParam('xDeviceToPixels');
    Inicoord.y.pix.onap=Inicoord.y.angle.onap * mglGetParam('yDeviceToPixels');
    
    %Calculate coordinates on aperture relative to screen's root.
    Inicoord2root.x.pix.onap=Inicoord.x.pix.onap + fixStimulus.pospix(1);
    Inicoord2root.y.pix.onap=Inicoord.y.pix.onap + fixStimulus.pospix(2);
    
    %Position mouse
    mglSetMousePosition(Inicoord2root.x.pix.onap, Inicoord2root.y.pix.onap, myscreen.screenNumber);
    
    %get Position (for when motion with 1/2 key presses)
    myscreen.mi = mglGetMouse(myscreen.screenNumber); %(pixels)
end

%confirmation
if (task.thistrial.thisseg == 4)
    
    %case subject move the response arrow with the powermate wheel, the
    %mouse or the cni touchpad
    if myscreen.responseType ~= 3
        %if no response recorded
        if task.thistrial.gotResponse ==0
            task.thistrial.prodcoor = [NaN NaN];
            task.thistrial.proddeg  = NaN;
            task.thistrial.reactionTime  = NaN;
            %if responses
        else
            %Get mouse position
            mi = mglGetMouse(myscreen.screenNumber); %get pixel positions.
            mouseinfo.x.pix = mi.x;
            mouseinfo.y.pix = mi.y;
            
            %Check that subject confirmed his choice: (keyboard press "1" down or mouse click")
            %if you want to you use mouse click to enter choice
            if myscreen.RespConfirm == 2
                mouseinfo.buttons = mi.buttons;
            elseif myscreen.RespConfirm == 1
                %case you use keyboard button '1' to confirm choice
                mouseinfo.buttons = mglGetKeys(19);
            end
            
            %Position mouse on aperture by moving mouse horizontally
            [x y] = Pos2apByMovingHoriz(mouseinfo.x.pix,task);
            mouseinfo.y.angle.onap = y.angle.onap;
            mouseinfo.x.angle.onap = x.angle.onap;
            
            %Position mouse on aperture
            %mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo);
            
            %Interface PowerMate
            %mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
            
            %Collect response
            task.thistrial.prodcoor = [mouseinfo.x.angle.onap mouseinfo.y.angle.onap];%(visual angle)
            [~,task.thistrial.proddeg] = SLcart2polar(task.thistrial.prodcoor);
            
        end
        %case response arrow is moved with keypresses "1"/"2"
        %and the reported estimate  is registered when confirmation key "3"
        %is pressed
    elseif myscreen.responseType == 3
        %if no response recorded
        if task.thistrial.gotResponse ==0
            task.thistrial.prodcoor = [NaN NaN];
            task.thistrial.proddeg  = NaN;
            task.thistrial.reactionTime  = NaN;
            %if response
        else
            %Position mouse on aperture
            %mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo);
            %Interface PowerMate
            %mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
            %Collect response
            task.thistrial.prodcoor = [myscreen.mouseinfo.x.angle.onap myscreen.mouseinfo.y.angle.onap];%(visual angle)
            [~,task.thistrial.proddeg] = SLcart2polar(task.thistrial.prodcoor);
        end
    end
    
    %print task info
    fprintf('-----------------------------------------------\n')
    fprintf('%s  %s  %s  %s  %s \n','trial','coh','direction','estimate','RT')
    fprintf('-----------------------------------------------\n')
    fprintf('  %i    %.2f    %i       %i     %.2f \n \n',...
        task.trialnum,...
        task.thistrial.myRandomCoh,...
        task.thistrial.myRandomDir,...
        round(task.thistrial.proddeg),...
        round(task.thistrial.reactionTime*100)/100)
end

%function that gets called to draw the stimulus each frame
function [task,myscreen] = updateScreenCallback(task, myscreen)

global stimulus;
global fixStimulus;
mglClearScreen;

%motion (outer annulus within which dots move while inner part
%with no dots and background color (black))
if (task.thistrial.thisseg == 2)
    stimulus = updateDots(stimulus,myscreen);
    %foveal black
    mglGluDisk(0,0,stimulus.rmaxResp,0,24,2)  %foveal black disk
    mglGluDisk(0,0,0.2,0.5,24,2)              %fixation dot
end
%response
if (task.thistrial.thisseg == 3)
    
    %draw resp arrow only foveal (separated from peripheral motion)
    mglGluDisk(0,0,2,0,24,2)                  %foveal black disk (deg)
    mglGluDisk(0,0,0.2,0.5,24,2)              %fixation (deg)
    
    %update until response
    if myscreen.responseType ~= 3
        if task.thistrial.gotResponse == 0
            
            %Position mouse on aperture
            mi = mglGetMouse(myscreen.screenNumber); %(pixels)
            mouseinfo.x.pix = mi.x;
            mouseinfo.y.pix = mi.y;
            
            %Position mouse on aperture by moving mouse horizontally
            [x y] = Pos2apByMovingHoriz(mouseinfo.x.pix,task);
            mouseinfo.y.angle.onap = y.angle.onap;
            mouseinfo.x.angle.onap = x.angle.onap;
            
            %Position mouse on aperture
            %mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo);
            
            %Interface PowerMate
            %mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
            
            %draw a white arrow (radius)
            mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
                mouseinfo.x.angle.onap,...
                mouseinfo.y.angle.onap,1,[1 1 1],1);
        end
        %case move with 1/2 keys and confirm with 3
        %once confirmed screen turned black until the end of the response
        %period
    elseif myscreen.responseType == 3 && task.thistrial.buttonConfstate==0
        
        %Position mouse on aperture
        mouseinfo.x.pix = myscreen.mi.x;
        mouseinfo.y.pix = myscreen.mi.y;
        
        %backup for motion
        myscreen.mouseinfo = mouseinfo;
        
        %once we start collecting response move
        %arrow cw or ccw
        if task.thistrial.buttonState(1)==1
            myscreen.mi.x = myscreen.mi.x + 1;
            %move cw ("2")
        elseif task.thistrial.buttonState(2)==1
            myscreen.mi.x = myscreen.mi.x - 1;
        end
        mouseinfo.x.pix = myscreen.mi.x;
        mouseinfo.y.pix = myscreen.mi.y;
        %Position mouse on apertut()re by moving mouse horizontally
        [x y] = Pos2apByMovingHoriz(mouseinfo.x.pix,task);
        myscreen.mouseinfo.y.angle.onap = y.angle.onap;
        myscreen.mouseinfo.x.angle.onap = x.angle.onap;
        %re-set segment response status to 0
        task.thistrial.gotResponse = 0;
        %draw a white arrow (radius)
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            myscreen.mouseinfo.x.angle.onap,...
            myscreen.mouseinfo.y.angle.onap,1,[1 1 1],1);
        
        %check whether response was confirmed in previous
        %frame and register response
        if task.thistrial.buttonState(3)==1
            task.thistrial.buttonConfstate = 1;
            task.thistrial.gotResponse = 1;
            %jump to confirmation at response
            if task.thistrial.gotResponse == 1
                task = jumpSegment(task,4);
            end
        end
    end
end

%confirmation
if (task.thistrial.thisseg == 4) %confirmation (segment 4)
    %case response arrow was moved with powermate wheel, mouse or cni
    %touchpad
    if myscreen.responseType ~= 3
        %translate mouse 2 aperture
        mi = mglGetMouse(myscreen.screenNumber); %deliver pixel positions, not angle
        mouseinfo.x.pix = mi.x;
        mouseinfo.y.pix = mi.y;
        
        %if you want to you use mouse click to enter choice
        mouseinfo.buttons = mi.buttons;%mglGetKeys(50); %mi.buttons;
        
        %Position mouse on aperture
        mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
        %Interface PowerMate
        mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
        
        %draw resp arrow only foveal (separated from peripheral motion)
        mglGluDisk(0,0,2,0,24,2)                  %foveal black disk (deg)
        mglGluDisk(0,0,0.2,0.5,24,2)              %fixation (deg)
        
        %Confirm subjects' choice by drawing a red arrow
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            mouseinfo.x.angle.onap,...
            mouseinfo.y.angle.onap,1,[1 0 0],1);
        
        %case response arrow is moved with keys "1"/"2" and response is
        %confirmed with "3"
    elseif myscreen.responseType == 3
        %draw resp arrow only foveal (separated from peripheral motion)
        mglGluDisk(0,0,2,0,24,2)                  %foveal black disk (deg)
        mglGluDisk(0,0,0.2,0.5,24,2)              %fixation (deg)
        
        %Confirm subjects' choice by drawing a red arrow
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            myscreen.mouseinfo.x.angle.onap,...
            myscreen.mouseinfo.y.angle.onap,2,[1 0 0],1);
    end
end

%feedback (foveal green arrow )
if (task.thistrial.thisseg == 5)
    
    %arrow coordinates
    coord = polar2cartesian(task.thistrial.myRandomDir,stimulus.rmaxResp);
    feedb.x.angle.onap = coord.x;
    feedb.y.angle.onap = coord.y;
    
    %draw resp arrow only foveal (separated from peripheral motion)
    mglGluDisk(0,0,2,0,24,2)    %foveal black disk (deg)
    mglGluDisk(0,0,0.2,0.5,24,2)%fixation (deg)
    
    %green feedback arrow
    mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
        feedb.x.angle.onap,...
        feedb.y.angle.onap,2,[0 1 0],1);
end

%function that gets called to get response
function [task,myscreen] = getResponseCallBack(task,myscreen)

%function to init the dot stimulus
function stimulus = initDots(stimulus,myscreen)

%convert the passed in parameters to real units
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax=15;end %(radius in angle)
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
function stimulus = updateDots(stimulus,myscreen)

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
function [mouseinfo] = pos2ap(myscreen,fixStimulus,mouseinfo,stimulus)

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
function [mouseinfo] = PowerMate2ap(mouseinfo,stimulus,task)
global Inicoord;
%Set speed, initial position and size of arrow.
wheelspeed=0.008;
%fixed
%initAngle=(225+270)/2*pi/180;
%random
initAngle.rad=task.thistrial.initAngledeg*pi/180;

r=stimulus.dots.rmax;
%Calculate arrow's coordinates on aperture
mouseinfo.x.angle.onap = r * cos(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
mouseinfo.y.angle.onap = r * sin(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap));
%function to convert from PowerMate's to aperture's coordinates
function [x y] = Pos2apByMovingHoriz(xpix,task)
global Inicoord;
global stimulus;
speed = stimulus.respSpeed;%Set speed,
initAngle.rad = task.thistrial.initAngledeg*pi/180;%random position
r = stimulus.rmaxResp;
%Calculate arrow's coordinates on aperture and move the mouse horizontally
%(we only operate on x coordinate) to adjust the response arrow position.
y.angle.onap = r*sin(initAngle.rad+speed*(xpix - Inicoord.x.pix.onap));
x.angle.onap = r*cos(initAngle.rad+speed*(xpix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
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
%fixationTrain
function [task myscreen] = fixationTrain(myscreen)

%fixationTrain.m
%
%       $Id: fixStairInitTask.m 708 2010-03-01 03:38:19Z justin $
%     usage: [fixTask myscreen]=fixStairInitTask(myscreen)
%        by: justin gardner
%      date: 09/07/06
% copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%   purpose: Implements a fixationTrain task. In this task, the fixationTrain cross
%            starts out cyan and then darkens twice. The fixationTrain cross then
%            turns yellow to indicate the response interval, and the subject
%            is required to press 1 or 2 to indicate in which interval the cross
%            appeared darker. The cross will then turn green or red to indicate
%            correct or incorrect responses. The dimness of the target is 
%            controlled by a 2 down 1 up staircase to keep task difficulty
%            the same. 
%
%            See testExperiment.m for how this is used in a task. If you want
%            to change parameters, before you call fixStairInitTask, set 
%            appropriate fields of the global variable fixStimulus. e.g.:
%
%            global fixStimulus
%            fixStimulus.interTime=1;
%
%            See the code, for a list of all parameters that can be changed.
%Modified by A Meso May 2012 for a fixationTrain synched to stimuli but without
%a related task to solve

%Modified by Steeve Laquitaine August 2012 for a fixationTrain synched to
%stimuli.

%Notes:
<<<<<<< HEAD
%Fixation point is a disk of 0.2 degrees radius in visual angle.

=======
%Fixation point is a disk of 0.2 degrees radius in visual angle.3111233
>>>>>>> affe4605afaf6eebdf2ef772d10ff2b1ff5d5601

%check arguments
if ~any(nargin == [1])
  help fixationTrain
  return
end

%create the stimulus for the experiment, use defaults if they are
%not already set
global fixStimulus;
myscreen = initStimulus('fixStimulus',myscreen);

%fixationTrain parameters
fixStimulus.fixWidth = 3; %(cm) at this size, fixationTrain also acts as a cardinal axis reference... 
fixStimulus.fixLineWidth = 100; %degree
fixStimulus.pos = [0 0]; %positon in the center of the screen

%init the task
%Before the 131220 this seglen was set at inf; but it nows produces an
%error. It may be due to mgl functions update. I set it now at a large
%number it seems to work fine.
task{1}.seglen = 100000;
[task{1} myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback);
<<<<<<< HEAD
=======

%function that gets called at the start of each segment
function [task myscreen] = fixStartSegmentCallback(task, myscreen)

global fixStimulus;
%if this is the inter stimulus interval
  %training mode, clear screen here
%  if fixStimulus.trainingMode,mglClearScreen;end 
  fixStimulus.thisColor=127; %grey so that it is not too strong over many trials
%function that gets called every frame udpate to draw the fixationTrain 
function [task myscreen] = fixDrawStimulusCallback(task, myscreen)
global fixStimulus;
>>>>>>> affe4605afaf6eebdf2ef772d10ff2b1ff5d5601

%Stimulus size (degree of visual angle)
mglGluDisk(0,0,0.2,0.5,24,2)
