% sltaskDual2AFC2.m
%
%
%      usage: sltaskDual2AFC
%         by: steeve laquitaine
%       date: 09/12/2015
%    purpose: Baruni et al., 2015 Nat Neuro.
%
%usage:
%
%   sltaskDual2AFC('displayName=Test')
%
%
%Description:
%
%- ~3 sec trial
%- 0.3s fixation (Baruni)
%- Two 0.5 deg radius reward cues at 2.25 deg from fixation (Baruni)
%- 2.5 deg radius, 2 cpd gratings at 6 deg from fixation (Pestilli)
%- 0.28 deg response cues (white line)
%
%To do 
%
% - add reward feedback: sounds
% - add masks


function myscreen = sltaskDual2AFC(varargin)

%check arguments
if ~any(nargin == [0 1])
    help sltaskDual2AFC
    return
end

%get inputs Default response device is powermate
displayName=[];
getArgs(varargin,{'displayName'});


%---------- init screen ----------
%displayName
myscreen.displayname = displayName;
%init
myscreen = initScreen(myscreen);
%hide cursor
mglDisplayCursor(0);

%init the stimulus
clear global stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

%these are the reserved colors, if you need them later
%you can display them by setting your color to the appropriate
%index in stimulus.colors.reservedColor e.g. to get the
%second color, in this case white, you would do
%mglClearScreen(stimulus.colors.reservedColor(2));
stimulus.colors.reservedColors = [0 0 0; 1 1 1; 0 0 1; 1 0.5 0];

%grating parameters
stimulus.grating.sf   =  2; %grating sf (cpd, Pestilli)
stimulus.grating.ecc  =  6; %ecc. in deg (Pestilli)
stimulus.grating.xpos = -stimulus.grating.ecc/sqrt(2); 
stimulus.grating.ypos = -stimulus.grating.ecc/sqrt(2);
stimulus.grating.windowType = 'gabor'; %should be gabor or thresh
stimulus.grating.width  = 5; %gaus mask (deg, Pestilli)
stimulus.grating.height = 5; %gaus mask (deg, Pestilli)
stimulus.grating.sdx = stimulus.grating.width/7;%gauss mask x std
stimulus.grating.sdy = stimulus.grating.width/7;%gauss mask y std
 

%-------------------- init reward cues -------------
%Two disks top right / top left
stimulus.rcue.rdius = 0.5; %(deg, Baruni)

%eccent. 2.25 deg eccent.
%note: ecc is the diagonal of a square
%with sides xpos and ypos
%WARNING ! this only works when xpos = ypos !!
stimulus.rcue.ecc = 2.25; %Baruni
stimulus.rcue.xpos = -stimulus.rcue.ecc/sqrt(2); %eccentricity in vis. ang (deg)
stimulus.rcue.ypos = -stimulus.rcue.ecc/sqrt(2);
stimulus = initrcues(stimulus);%init Cues

%---------------------init task -----------------
%trial segments fix - rew cues - delay/stream - target/dist - mask - resp
%- print out phase 
%no grating stream for now
%Baruni
%task{1}{1}.segmin = [.3 .25 .35+.5  .1 .06 3 0.01];
%task{1}{1}.segmax = [.3 .25 .35+1.5 .15 .06 3 0.01];
task{1}{1}.segmin = [.3 .25 .35+.5  .1 .06 3 0.01];
task{1}{1}.segmax = [.3 .25 .35+1.5 .15 .06 3 0.01];


%set up task
%Orientations range from 0 to 90 deg
task{1}{1}.randVars.orientG1  = 0:3:90; %bot grating
task{1}{1}.randVars.orientG2  = 90:-3:0;%top grating
task{1}{1}.parameter.contrast = 1;
task{1}{1}.numTrials = length(task{1}{1}.randVars.orientG1);
rndT = randperm(task{1}{1}.numTrials);
rndT2 = randperm(task{1}{1}.numTrials);

%randomize G1,G2 orientations
task{1}{1}.randVars.orientG1  = task{1}{1}.randVars.orientG1(rndT);
task{1}{1}.randVars.orientG2  = task{1}{1}.randVars.orientG2(rndT2);

%get response and reaction time
%response segment is set at "3" or "4" for
%mouse event and 1 for keyboard press.
%We use keyboard press "1" and "2"
%for "near horisontal"/"vertical" responses
task{1}{1}.getResponse = [0 0 0 0 0 1];

%init response
task{1}{1}.randVars.calculated.myresp = nan;

%----------------- init stimulus ---------------
stimulus = initGratings(stimulus,myscreen,task);

%-------------- init response cues -------------
nTrials = task{1}{1}.numTrials;
task{1}{1}.randVars.respCue = (rand(nTrials,1)>0.5)+1;

%init tasks
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);
end

%run the eye calibration
myscreen = eyeCalibDisp(myscreen);

%run the tasks
tnum = 1;
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

%if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%reward cues
function stimulus = initrcues(stimulus)

%set max color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

%set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
for i = 1:stimulus.colors.nReservedColors
    stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
end

%get gray luminance
stimulus.colors.nGratingColors = maxIndex+1-stimulus.colors.nReservedColors;
stimulus.colors.minGratingIndex = maxIndex+1 - stimulus.colors.nGratingColors;
stimulus.colors.midGratingIndex = stimulus.colors.minGratingIndex+floor(stimulus.colors.nGratingColors/2);
stimulus.colors.grayColor = stimulus.colors.midGratingIndex/maxIndex;

%a simple window
disk1 = slmakeDisk(stimulus.rcue.rdius,0,0);

%rcue top-right : make indices from 0 to 255
rcue1 = nan(size(disk1));
rcue1(disk1==0) = stimulus.colors.midGratingIndex;
rcue1(disk1==1) = 3;
stimulus.rcue1 = mglCreateTexture(rcue1);

%rcue bottom-left : make indices from 0 to 255
rcue2 = nan(size(disk1));
rcue2(disk1==0) = stimulus.colors.midGratingIndex;
rcue2(disk1==1) = 2;
stimulus.rcue2 = mglCreateTexture(rcue2);



%Gratings
function stimulus = initGratings(stimulus,myscreen,task)

%set max color index (for 24 bit color we have 8 bits per channel, so 255)
maxIndex = 255;

%get gamma table
if ~isfield(myscreen,'gammaTable')
    stimulus.linearizedGammaTable = mglGetGammaTable;
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
    disp(sprintf('(taskTemplateContrast10bit:initGratings) No gamma table found in myscreen. Contrast'));
    disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
    disp(sprintf('         for this monitor.'));
    disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

disppercent(-inf,'Creating grating textures');

%calculate colors info
%# of reserved colors
stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);

%# of colors possible for gratings, make sure that we
%have an odd #
stimulus.colors.nGratingColors = maxIndex+1-stimulus.colors.nReservedColors;
if iseven(stimulus.colors.nGratingColors)
    stimulus.colors.nGratingColors = stimulus.colors.nGratingColors-1;
end
%min, mid and max index of gratings colors (index values are 0 based)
stimulus.colors.minGratingIndex = maxIndex+1 - stimulus.colors.nGratingColors;
stimulus.colors.midGratingIndex = stimulus.colors.minGratingIndex+floor(stimulus.colors.nGratingColors/2);
stimulus.colors.maxGratingIndex = maxIndex;
%# of contrasts we can display (not including 0 contrast)
stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nGratingColors/2);

%get the color value for gray (i.e. the number between 0 and 1 that corresponds to the midGratingIndex)
stimulus.colors.grayColor = stimulus.colors.midGratingIndex/maxIndex;

%set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
    stimulus.colors.reservedColor(i) = (i-1)/maxIndex;
end

%make the window through with the gratings will be displayed
gaussianWin = mglMakeGaussian(stimulus.grating.width,stimulus.grating.height,stimulus.grating.sdx,stimulus.grating.sdy);
if strcmp(stimulus.grating.windowType,'gabor')
    %a gaussian window
    win = maxIndex-maxIndex*gaussianWin;
else
    %a simple window
    win = maxIndex-maxIndex*(gaussianWin>exp(-1/2));
end
mask = ones(size(win,1),size(win,2),4)*stimulus.colors.midGratingIndex;
mask(:,:,4) = win;
stimulus.mask = mglCreateTexture(mask);

%Always max contrast
iContrast = stimulus.colors.nDisplayContrasts;

%params
Gwidth = stimulus.grating.width;
Gheight = stimulus.grating.height;
Gsf = stimulus.grating.sf;
GmidIx = stimulus.colors.midGratingIndex;

%all gratings
%bottom-left gratings
angG1 = task{1}{1}.randVars.orientG1;
nAng = length(angG1);
for i = 1 : nAng
    disppercent(i/stimulus.colors.nDisplayContrasts);
    if myscreen.userHitEsc,mglClose;keyboard,end
    %make grating
    thisGrating1 = round(iContrast*mglMakeGrating(Gwidth,Gheight,Gsf,angG1(i),0)+GmidIx);
    %create texture
    stimulus.tex1(i) = mglCreateTexture(thisGrating1);
end

%top-right gratings
angG2 = task{1}{1}.randVars.orientG2;
nAng = length(angG2);
for i = 1 : nAng
    disppercent(i/stimulus.colors.nDisplayContrasts);
    if myscreen.userHitEsc,mglClose;keyboard,end
    %make grating
    thisGrating2 = round(iContrast*mglMakeGrating(Gwidth,Gheight,Gsf,angG2(i),0)+GmidIx);
    %create texture
    stimulus.tex2(i) = mglCreateTexture(thisGrating2);
end


%----------------------- each segment ----------------------
function [task myscreen] = startSegmentCallback(task,myscreen)

%set the max contrast we can display
%disp(sprintf('(taskTemplateContrast10bit:startSegmentCallback) Displaying contrast of %f',task.thistrial.contrast));
setGammaTableForMaxContrast(task.thistrial.contrast);

%print out
if task.thistrial.thisseg == 7
    if ~any(task.thistrial.buttonState)
        %missing responses
        if ~any(task.thistrial.myresp)==1
            task.thistrial.myresp = nan;
        end
        %get reaction time
        if ~isfield(task.thistrial,'reactionTime')
            task.thistrial.reactionTime = NaN;
        end
    end
    %print trial # with key
    fprintf(' %i    %i   %i   %i  %.1f \n',task.trialnum,...
        task.thistrial.orientG1,...
        task.thistrial.orientG2,...
        task.thistrial.myresp,...
        task.thistrial.reactionTime)
end

%----------------------- each frame -------------------
function [task myscreen] = updateScreenCallback(task,myscreen)

global stimulus;

%clear screen to gray
mglClearScreen(stimulus.colors.grayColor);

%get the contrast index
%contrastIndex = getContrastIndex(task.thistrial.contrast);

%blt texture
%position
xpos = stimulus.grating.xpos;
ypos = stimulus.grating.ypos;

%rcue pos
xCue = stimulus.rcue.xpos;
yCue = stimulus.rcue.ypos;

%reward cue phase
if task.thistrial.thisseg == 2
    %cue1: top right
    mglBltTexture(stimulus.rcue1,[-xCue -yCue]);
    %cue 2: bottom left
    mglBltTexture(stimulus.rcue2,[xCue yCue]);
end

%Target/Distractor phase
if task.thistrial.thisseg == 4
    %draw grating 1
    mglBltTexture(stimulus.tex1(task.trialnum),[xpos ypos stimulus.grating.height]);
    %gaussian mask 1
    mglBltTexture(stimulus.mask,[xpos ypos]);
    
    %draw grating 2
    mglBltTexture(stimulus.tex2(task.trialnum),[-xpos -ypos stimulus.grating.height]);
    %gaussian mask 2
    mglBltTexture(stimulus.mask,[-xpos -ypos]);
end

%response cue
if task.thistrial.thisseg == 6    
    %top right grating
    if task.thistrial.respCue == 1    
        mglLines2(.2,.2,.4,.4,1,stimulus.colors.reservedColor(2),0);
        %bottom left grating
    elseif task.thistrial.respCue == 2
        mglLines2(-.2,-.2,-.4,-.4,1,stimulus.colors.reservedColor(2),0);
    end
end

%put up the fixation cross with reserved color 5
mglFixationCross(0.5,1,stimulus.colors.reservedColor(2));

%---------------------- get response ------------------------
function [task,myscreen] = getResponseCallBack(task,myscreen)

%response
response = find(task.thistrial.buttonState);
task.thistrial.myresp = response;

%When subject chose, jump
task = jumpSegment(task);



%control gamma contrast
function setGammaTableForMaxContrast(maxContrast)

global stimulus;
% if you just want to show gray, that's ok, but to make the
% code work properly we act as if you want to display a range of contrasts
if maxContrast <= 0,maxContrast = 0.01;end

% set the reserved colors
gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the gamma table
if maxContrast > 0
    % create the rest of the gamma table
    cmax = 0.5+maxContrast/2;cmin = 0.5-maxContrast/2;
    luminanceVals = cmin:((cmax-cmin)/(stimulus.colors.nGratingColors-1)):cmax;
    
    % now get the linearized range
    redLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
    
    % add these values to the table
    gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,:)=[redLinearized;greenLinearized;blueLinearized]';
else
    % if we are asked for 0 contrast then simply set all the values to gray
    gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,1)=interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,0.5,'linear');
    gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,2)=interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,0.5,'linear');
    gammaTable((stimulus.colors.minGratingIndex:stimulus.colors.maxGratingIndex)+1,3)=interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,0.5,'linear');
end

% set the gamma table
mglSetGammaTable(gammaTable);

% remember what the current maximum contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

