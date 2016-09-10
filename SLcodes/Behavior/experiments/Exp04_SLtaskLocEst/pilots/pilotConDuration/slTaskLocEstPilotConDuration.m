%slTaskLocEstPilotConDuration.m
%
%
% author: Steeve Laquitaine
%   date: 151109
%purpose: location estimation task
% status: completed (fine contrast)
%
%  usage:
%
%           slTaskLocEst('fMRI')
%           slTaskLocEst('psychophysics','displayName=Test','responseObject=dot','responseDevice=1')
%           slTaskLocEst('fMRI','params=steeve_exp12_metho_Pstd010_mean225_con006012024_dir36_t107_074_033perCoh_controlRndInitPosSymPrior_131224')
%           slTaskLocEst('psychophysics','testOneContrast',0.3)
%              
%Inputs
%
%                 'psychophysics': ('fMRI') psychophysics/fmri task
%
%                'responseDevice': 'responseDevice=0': (mouse, 1:powermate, 2: keyboard(not working yet))
%                   'displayName': 'displayName=Test'
%                'responseObject': 'responseObject=dot'('arrow') displays a dot/arrow subject can move around to report wedge location
%       
%
%
%ref: 
%   taskTemplateContrast10bit.m
%   http://gru.stanford.edu/doku.php/mgl/taskReferenceHowTos#how_to_use_10-bit_contrast
%   http://gru.stanford.edu/doku.php/mgl/functionreferencegammatable

function myscreen = slTaskLocEstPilotConDuration(varargin)

%get args
eval(evalargs(varargin,0,0,{'displayName'}));

%-----------------------------init Task----------------------------
%Main task
%Psychophysics
%same as dot direction task
if ~ieNotDefined('psychophysics')
    %trial
    task{2}{1}.segmin = [1 .3 5 .1 .1];
    task{2}{1}.segmax = [1 .3 5 .1 .1];
    myscreen.myExperiment = 'psychophysics';    
    fprintf('(slTaskLocEst) Psychophysics task. \n')    
    %default
    if ieNotDefined('params')    
        fprintf(['\n \n ------ ! WARNING !------ You did not enter task ' ...
                 'params. Task is running in INSTRUCTION MODE !!. \n \n'])
        YoN = input(['(slTaskLocEst) No task params. Should I use ' ...
                     'defaultparams ? (y/n) :'],'s');    
        if strcmp(YoN,'y')
            params.task.parameter.loc.series = repmat([5:10:355]',3,1);
            numtrials = numel(params.task.parameter.loc.series); %location
            params.task.parameter.loc.trial     = 1:numtrials;    
            params.task.parameter.loc.con = repmat([1;0.5;0.2],36,1);%contrast
            params.task.parameter.loc.mean  = 225;
            params.task.parameter.loc.modes = 225;
            params.task.parameter.loc.std   = inf;
            params.task.parameter.loc.strength = 0;
            params.task.parameter.loc.sample.degree = unique(params.task.parameter.loc.series);
        else
            keyboard
        end
    else
        %load params
        if ~ieNotDefined('params')        
            myPath = mfilename('fullpath');
            cd(fileparts(myPath))
            %Load parameters
            params = load(['TaskParameters/',params,'.mat']);
        end
    end    
elseif ~ieNotDefined('fMRI')
    %case fMRI task
    %1s fixation + stim + 11s fixation ITI
    task{2}{1}.segmin = [1 .3 5 .1 .1 8.8];
    task{2}{1}.segmax = [1 .3 5 .1 .1 8.8];
    myscreen.myExperiment = 'fMRI';
    fprintf('fMRI task. \n')
    %default
    if ieNotDefined('params')    
        fprintf(['\n \n ------ ! WARNING !------ You did not enter task ' ...
                 'params. Task is running in INSTRUCTION MODE !!. \n \n'])
        YoN = input(['(slTaskLocEst) No task params. Should I use default ' ...
                     'params ? (y/n)','s']);    
        if strcmp(YoN,'y')
            params.task.parameter.loc.series = repmat([5:10:355]',3,1);
            numtrials = numel(params.task.parameter.loc.series); %location
            params.task.parameter.loc.trial     = 1:numtrials;    
            params.task.parameter.loc.con = repmat([1;0.5;0.2],36,1);%contrast
            params.task.parameter.loc.mean  = 225;
            params.task.parameter.loc.modes = 225;
            params.task.parameter.loc.std   = inf;
            params.task.parameter.loc.strength = 0;
            params.task.parameter.loc.sample.degree = unique(params.task.parameter.loc.series);
        else
            keyboard
        end
    else
        %load params
        if ~ieNotDefined('params')        
            myPath = mfilename('fullpath');
            cd(fileparts(myPath))
            %Load parameters
            params = load(['TaskParameters/',params,'.mat']);
        end
    end
else
    fprintf(['Is the task "fMRI" or "psychophysics" ? Please rerun ' ...
             'with something like slTaskLocEst("fMRI") or "psychophysics" \n'])
    mglClose
    keyboard
end

%-------------------- Screen -----------------------------

%mandatory inputs
if ieNotDefined('displayName'),
    displayName = input('input displayName (e.g., test): ','s');
end
if ieNotDefined('responseObject'),
    responseObject = input(['input responseObject ("arrow" or "dot") ' ...
                        ':'],'s');
end
myscreen.responseObject = responseObject;
if ieNotDefined('responseDevice'),
    myscreen.responseDevice = input(['responseDevice (0: mouse, 1: ' ...
                        'powermate, keyboard:2) :']);
else 
    myscreen.responseDevice = responseDevice;
end
    
%save data in user dir data/SLtaskLocEst
myscreen.saveData = 1;
myscreen.displayname = displayName;
myscreen = initScreen(myscreen);

%----------------------- stimulus -----------------------
global stimulus;

%init gamma table
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%Colors & contrasts
stimulus.colors.resvdTriplet = [0 0 0; 1 1 1; 0 1 0; 1 0 0];
stimulus.colors.nresvd = size(stimulus.colors.resvdTriplet,1);

%# of possible stim. colors (e.g., 255 - nres)
maxIx = 255;
stimulus.colors.nPossCol = maxIx+1 - stimulus.colors.nresvd;
if iseven(stimulus.colors.nPossCol)
  stimulus.colors.nPossCol = stimulus.colors.nPossCol - 1;
end

%min, mid and max index row index of stim colors in gammaTable 
%(index values are 0 based). e.g., 
%min/mid/max/StimIx = 253/254/255;
stimulus.colors.minStimIx = maxIx+1 - stimulus.colors.nPossCol;
stimulus.colors.midStimIx = stimulus.colors.minStimIx + floor(stimulus.colors.nPossCol/2);
stimulus.colors.mxSIx     = maxIx;
stimulus.colors.grayColor = stimulus.colors.midStimIx/maxIx;
stimulus.background       = stimulus.colors.grayColor;

%set the reserved colors as normalized 
%luminance indices (always [0 [1:255]/255))
%Those indices are associated with each of the
%255 rows in the gammaTable and with each of the
%corresponding luminances triplets for RGB.
%e.g., 0 is always associated with the gammaTable'
%first RGB luminance triplet (1st row) and 1 with the last 
%(255th row)%
% e.g., the indices
%of 5 reserved colors positioned at the first rows
%of the gamma table will be [0:4]/255
%Those are the indices called to display their associated 
%luminances.
for i = 1 : stimulus.colors.nresvd
  stimulus.colors.resvdNormIx(i) = (i-1)/maxIx;
end

%#of contrasts we can display (not including 0 contrast)
stimulus.colors.nDisplayContrasts = floor(stimulus.colors.nPossCol/2);

%settings that are used to adjust the position on the screen
%the stimuli are shown in - for cases when the subject can
%not see the whole screen
%sizeStim: 1/36 allows 36 non overlapping location.
if ieNotDefined('imageWidth'),imageWidth = [];end
if ieNotDefined('imageHeight'),imageHeight = [];end
if ieNotDefined('sizeStim'),sizeStim = 1/36;end
if ieNotDefined('elementSize'),elementSize = [];end
if ieNotDefined('xOffset'),xOffset = 0;end
if ieNotDefined('yOffset'),yOffset = 0;end
if ~isfield(stimulus,'dots') || ~isfield(stimulus.dots,'rmax'), stimulus.dots.rmax=14;end %(radius in angle, max on VPIxx)

%screen size settings as default
if isempty(imageWidth)
    stimulus.imageWidth = myscreen.imageWidth;
else
    stimulus.imageWidth = imageWidth;
end
if isempty(imageHeight)
    stimulus.imageHeight = myscreen.imageHeight;
else
    stimulus.imageHeight = imageHeight;
end

%Wedge min/max radius size
stimulus.minRadius = 0;
stimulus.maxRadius = stimulus.dots.rmax;
stimulus.sizeStim  = sizeStim;

%Wedge elements
stimulus.elementAngleSize      = 5; %size in degrees
stimulus.elementRadiusSize     = 2; %radial length
stimulus.elementRadialVelocity = 0; %radial speed (0 for static)
stimulus.elementVelocity       = 6; %element velocity
stimulus.FlickerFreq           = 8;%8 %flickering frequency in hertz

%wedge flickering frequency. 
%#of frames during which stimulus is static
frameRate                    = mglGetParam('frameRate');
frameDuration                = 1/frameRate;
stimulus.FlickerPeriod       = 1/stimulus.FlickerFreq;
stimulus.nFrameStill         = round(frameRate*stimulus.FlickerPeriod);
stimulus.nFrameStillDuration = stimulus.nFrameStill * frameDuration;

%screen size settings
stimulus.xOffset = xOffset;
stimulus.yOffset = yOffset;

%init stim
myscreen = initStimulus('stimulus',myscreen);

%------------------------- fixation --------------------------

myscreen.fixationColor = stimulus.colors.resvdNormIx(2);
myscreen.fixRadius     = 0.2;
myscreen.fixPosX       = 0;
myscreen.fixPosY       = 0;
mglDisplayCursor(0);

%!!ensure that screen is not flipped
%this is important !! it changes the value of the displayed location
if ~sum(myscreen.flipHV)==0
    fprintf(['!! WARNING !! The screen is flipped.',...
        'This might flip the stimulus displayed location. \n'])
    fprintf('If you must flip the screen, check that the displayed location on screen matches the output location in the terminal. \n')
    mglClose
    keyboard
end

%fixation 
global fixStimulus;
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2;
[task{1},myscreen] = fixationLoc(myscreen);

%----------------------- Task ----------------------------------

%randomize trials
myTrialRand = randperm(numel(params.task.parameter.loc.series));

%if instruction 
if ieNotDefined('params') 
    myTrialRand = 1 : numel(params.task.parameter.loc.series);
end

%Task variables
task{2}{1}.randVars.myRandomloc = params.task.parameter.loc.series(myTrialRand);
task{2}{1}.randVars.myRandomCon = ...
    params.task.parameter.loc.con(myTrialRand);
if ~ieNotDefined('psychophysics')
    params.task.parameter.loc.trialType = zeros(length(params.task.parameter.loc.series),1);
elseif ~ieNotDefined('fMRI')
    params.task.parameter.loc.trialType = ones(length(params.task.parameter.loc.series),1);
end
task{2}{1}.randVars.trialType   = params.task.parameter.loc.trialType(myTrialRand);
task{2}{1}.randVars.myMean      = params.task.parameter.loc.mean;
if isfield(params.task.parameter.loc,'modes')
    task{2}{1}.randVars.myModes = params.task.parameter.loc.modes;
end
task{2}{1}.randVars.myStrength  = params.task.parameter.loc.std;
task{2}{1}.randVars.myLocations = params.task.parameter.loc.sample;

%if test one contrast
if slIsInput(varargin,'testOneContrast')
    testCon = varargin{find(strcmp(varargin,'testOneContrast'))+1};
    task{2}{1}.randVars.myRandomCon = ones(length(task{2}{1}.randVars.myRandomloc),1)*testCon;
    slPrintfStr('SLtaskLocEst',['Testing contrast ' num2str(testCon)])
end

%randomize initial mouse location
task{2}{1}.randVars.initAngledeg = randi([0,359],[numel(task{2}{1}.randVars.myRandomloc),1]);

%get response and reaction time
%response segment is set at "3" or "4" for mouse event
task{2}{1}.getResponse = [0 0 1 0 0];
task{2}{1}.numTrials = numel(task{2}{1}.randVars.myRandomloc);

%init backup
task{2}{1}.randVars.calculated.prodcoor = [nan nan];
task{2}{1}.randVars.calculated.proddeg = nan;

%------------------------- set gamma table ---------------------------

%set the max contrast we can display
%It is set to the max contrast input
stimulus.maxCon = max(task{2}{1}.randVars.myRandomCon);
disp(sprintf('(slTaskLocEst) Displaying contrast of %f',stimulus.maxCon));
setGammaTableForMaxContrast(stimulus.maxCon);

%set each trial actual contrast
t = mglGetGammaTable; %new gammaT
gammaTableNew = [t.redTable' t.greenTable' t.blueTable'];
for i = 1 : length(task{2}{1}.randVars.myRandomloc)
    [actualCon(i),LumMaxIx(i),LumMinIx(i)] = getNormIdxContrast(task{2}{1}.randVars.myRandomCon(i),gammaTableNew);
end
task{2}{1}.randVars.actualCon = actualCon;
task{2}{1}.randVars.LumMax = LumMaxIx;
task{2}{1}.randVars.LumMin = LumMinIx;

%--------------------------- start ------------------------------------

%Init task
[task{2}{1},myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);

%Init stimulus
stimulus = initWedge(stimulus,myscreen,task);

%run eye calibration
myscreen = eyeCalibDisp(myscreen);

%Main display loop
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
    
    %update the dots
    [task{2},myscreen,phaseNum] = updateTask(task{2},myscreen,phaseNum);
    
    %update the fixation task
    [task{1},myscreen] = updateTask(task{1},myscreen,1);
    myscreen = tickScreen(myscreen,task);
end

%End
myscreen = endTask(myscreen,task);
clear global stimulus;
clear global fixStimulus;


%run at each trial segment
function [task,myscreen] = startSegmentCallback(task, myscreen)
global stimulus;
global fixStimulus;
global Inicoord;

%location
if (task.thistrial.thisseg == 2)   
    %update location strength each trial
    stimulus.dots.contrast = task.thistrial.actualCon;
    stimulus.dots.LumMax = task.thistrial.LumMax;
    stimulus.dots.LumMin = task.thistrial.LumMin;
    %actual contrast
    slPrintfStr('slTaskLocEst',[' Contrast(actual): ' num2str(stimulus.dots.contrast) '| (wanted) ' num2str(task.thistrial.myRandomCon)])
    %update mask location each trial
    stimulus.currentMask = stimulus.currentMask + 1;    
    %update Frame number each frame
    stimulus.countFrame = 0;    
    %warning when flickering period > stimulus duration
    if stimulus.nFrameStillDuration > task.thistrial.seglen(task.thistrial.thisseg)        
        minVisibleFlickerFreq = 1/(task.thistrial.seglen(task.thistrial.thisseg)/2);
        fprintf(['\n ------ ! WARNING ! ---- You flickering period is longer than the stimulus duration',...
            'You won t see flickering. Increase stimulus duration or flickering frequency. \n'])        
        fprintf([num2str(minVisibleFlickerFreq),'Hertz is the minimum visible flickering frequency',...
            ' with the current stimulus duration \n'])
        mglClose
        keyboard
    end
    stimulus.loc = task.thistrial.myRandomloc;    
end
stimulus.loc = task.thistrial.myRandomloc;


%fMRI (stim + long ITI)
if task.thistrial.trialType == 1
    
    %psyphy: stim+  "response", "confirmation" and "feedback"
elseif task.thistrial.trialType == 0
    
    %response
    if (task.thistrial.thisseg == 3)
        %Position mouse in the center of the screen
        %mglSetMousePosition(fixStimulus.pospix(1), fixStimulus.pospix(2), myscreen.screenNumber);
        %Calculate angle coordinates on aperture relative to center
        coord = polar2cartesian(task.thistrial.initAngledeg,stimulus.dots.rmax);
        Inicoord.x.angle.onap = coord.x;
        Inicoord.y.angle.onap = coord.y;
        %Calculate pixels coordinates on aperture relative to center
        Inicoord.x.pix.onap = Inicoord.x.angle.onap * mglGetParam('xDeviceToPixels');
        Inicoord.y.pix.onap = Inicoord.y.angle.onap * mglGetParam('yDeviceToPixels');
        %Calculate coordinates on aperture relative to screen's root.
        Inicoord2root.x.pix.onap = Inicoord.x.pix.onap + fixStimulus.pospix(1);
        Inicoord2root.y.pix.onap = Inicoord.y.pix.onap + fixStimulus.pospix(2);
        %Position mouse
        mglSetMousePosition(Inicoord2root.x.pix.onap, Inicoord2root.y.pix.onap, myscreen.screenNumber);
    end
    
    %confirmation
    if (task.thistrial.thisseg == 4)
        %Get mouse position
        mi = mglGetMouse(myscreen.screenNumber); %get pixel positions.
        mouseinfo.x.pix = mi.x;
        mouseinfo.y.pix = mi.y;
        %Check if subject confirmed his choice: ("space bar is down")
        %if you want to you use mouse click to enter choice
        %mouseinfo.buttons=mi.buttons;
        %if you want to you use keyboard button '1' to enter choice
        mouseinfo.buttons = mglGetKeys(19);
        %Position mouse on aperture
        mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
        %interface PowerMate
        if myscreen.responseDevice==1
            mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
        end
        %Collect response
        task.thistrial.prodcoor = [mouseinfo.x.angle.onap mouseinfo.y.angle.onap];%(visual angle)
        [~,task.thistrial.proddeg] = SLcart2polar(task.thistrial.prodcoor);
        %print trial number with the key and direction choosen.
        fprintf('%s  %s  %s  %s  %s \n','trial','trial type','contrast','polarAngle','response')
        fprintf('    %i    %i   %1.03f   %i   %1.01f \n \n',task.trialnum,...
            task.thistrial.trialType,...
            task.thistrial.myRandomCon,...
            task.thistrial.myRandomloc,...
            task.thistrial.proddeg)
        %then jump to feedback
        task = jumpSegment(task,4);
    end
    
    %feedback
    if (task.thistrial.thisseg == 5)
    end
    
    %ITI (fixation)
    if (task.thistrial.thisseg == 6)
    end
end

%run at each screen frame
function [task,myscreen] = updateScreenCallback(task, myscreen)
global stimulus;
global fixStimulus;
mglClearScreen(stimulus.colors.grayColor)

%need to be re-set after eyecalib
setGammaTableForMaxContrast(stimulus.maxCon);

%location
if (task.thistrial.thisseg == 2)    
    %update stimulus
    stimulus = updateWedge(stimulus,myscreen);    
    %update frame count each frame
    stimulus.countFrame = stimulus.countFrame + 1;    
    %update stimulus only when number of frames satisfying our flickering
    %frequency is reached otherwise stimulus remains frozen
    if SLisMultiple(stimulus.countFrame,stimulus.nFrameStill);
        stimulus.currentFlicker = stimulus.currentFlicker + 1;
    end
end

%case main trial, "ITI"
if task.thistrial.trialType == 1
    return    
    
elseif task.thistrial.trialType == 0
    %case catch trial, "response", "conf", "feedback"
    %response
    if (task.thistrial.thisseg == 3)        
        %Position mouse on aperture
        mi = mglGetMouse(myscreen.screenNumber); %(pixels)
        mouseinfo.x.pix=mi.x;
        mouseinfo.y.pix=mi.y;
        mouseinfo=pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
        %interface PowerMate
        if myscreen.responseDevice==1
            mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
        end
        %interface keyboard
        if myscreen.responseDevice==2          
        end
        %draw a white arrow (radius)
        if strcmp(myscreen.responseObject,'arrow');
            mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
                mouseinfo.x.angle.onap,...
                mouseinfo.y.angle.onap,2,stimulus.colors.resvdNormIx(2),1);
            %draw a white dot at the periphery (radius)
        elseif strcmp(myscreen.responseObject,'dot');
            mglPoints2(mouseinfo.x.angle.onap,mouseinfo.y.angle.onap,6,stimulus.colors.resvdNormIx(2))
        end
    end
    
    %confirmation
    if (task.thistrial.thisseg == 4)
        %Position mouse on aperture
        %deliver pixel positions, not angle
        mi = mglGetMouse(myscreen.screenNumber);
        mouseinfo.x.pix = mi.x;
        mouseinfo.y.pix = mi.y;
        mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo,stimulus);
        %Interface PowerMate
        if myscreen.responseDevice==1
            mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
        end
        %draw a red arrow (radius)
        if strcmp(myscreen.responseObject,'arrow');
            mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
                mouseinfo.x.angle.onap,...
                mouseinfo.y.angle.onap,2,stimulus.colors.resvdNormIx(4),1);
            %draw a red dot at the periphery (radius)
        elseif strcmp(myscreen.responseObject,'dot');
            mglPoints2(mouseinfo.x.angle.onap, ...
                       mouseinfo.y.angle.onap,10,stimulus.colors.resvdNormIx(4))
        end
    end
    
    %feedback
    if (task.thistrial.thisseg == 5)
        coord = polar2cartesian(task.thistrial.myRandomloc,stimulus.dots.rmax);
        feedbackcoord.x.angle.onap=coord.x;
        feedbackcoord.y.angle.onap=coord.y;
        %draw a black arrow (radius)
        if strcmp(myscreen.responseObject,'arrow');
            mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
                feedbackcoord.x.angle.onap,...
                feedbackcoord.y.angle.onap,4,stimulus.colors.resvdNormIx(1),1);
            %draw a black dot at the periphery (radius)
            %green is not visible on grey background
        elseif strcmp(myscreen.responseObject,'dot');
            mglPoints2(feedbackcoord.x.angle.onap,feedbackcoord.y.angle.onap,12,stimulus.colors.resvdNormIx(1))
        end
    end
    
    %ITI (fixation)
    if (task.thistrial.thisseg == 6)
    end
end

%get response
function [task,myscreen] = getResponseCallBack(task,myscreen)

%When subject chose, confirm (jump to segment 4 in updateScreenCallback)
task=jumpSegment(task,4);




%init stim
function stimulus = initWedge(stimulus,myscreen,task)

% round to nearest quarter of a degree, this reduces
% some edge effects
stimulus.maxRadius = floor(stimulus.maxRadius/.25)*.25;
disp(sprintf('Stimulus radius = [%0.2f %0.2f] degrees',stimulus.minRadius,stimulus.maxRadius));

% calculate some parameters
% size of wedges
stimulus.wedgeAngle = 360*stimulus.sizeStim;

% how much to step the wedge angle by
%stimulus.wedgeStepSize = 360/stimulus.stepsPerCycle;
sizeStim = stimulus.sizeStim;

%we only need to recompute the mglQuad points of the elements if something has
%changed in the stimulus. This is for the radial element pattern
if ~isfield(stimulus,'last') || ~isfield(stimulus,'x') || ...
        (stimulus.elementAngleSize ~= stimulus.last.elementAngleSize) || ...
        (stimulus.elementRadiusSize ~= stimulus.last.elementRadiusSize) || ...
        (stimulus.maxRadius ~= stimulus.last.maxRadius) || ...
        (stimulus.minRadius ~= stimulus.last.minRadius)
    
    %all the angles that the elements will be made up of
    allAngles = (0:stimulus.elementAngleSize:(360-stimulus.elementAngleSize));
    
    %all the phases. The phase refers to the radial position of the
    %black and white pattern (the pattern that is seen as moving
    %in the stimulus). There are two sets here since the wedges slide
    %against each other. That is every other sector will go in a
    %different direction.
    %case static stimulus
    if stimulus.elementRadialVelocity ==0
        allPhases1 = 1;
        allPhases2 = 1;
    else
        allPhases1 = 0:(stimulus.elementRadialVelocity/myscreen.framesPerSecond):(stimulus.elementRadiusSize*2);
        allPhases2 = fliplr(allPhases1);
    end
    
    %number of flickering
    allFlickers = task{2}{1}.segmax*myscreen.framesPerSecond;
    
    disppercent(-inf,'(mglRetinotopy) Calculating coordinates of elements in stimulus pattern');
%     for flickerNum = 1 : allFlickers
        for phaseNum = 1:length(allPhases1)
            stimulus.x{phaseNum} = [];stimulus.y{phaseNum} = [];stimulus.c{phaseNum} = [];
            
            %for each angle that the element contains
            for angleNum = 1:length(allAngles)
                
                %get the angle
                angle = allAngles(angleNum);
                
                % choose which phase we are going to be
                if isodd(angleNum)
                    thisMinRadius = stimulus.minRadius - allPhases1(phaseNum);
                else
                    thisMinRadius = stimulus.minRadius - allPhases2(phaseNum);
                end
                
                %all the radiuses
                allRadius = thisMinRadius : stimulus.elementRadiusSize : stimulus.maxRadius;
                
                %now create all the quads for this wedge
                for radiusNum = 1 : length(allRadius)
                    radius = allRadius(radiusNum);
                    
                    if (radius+stimulus.elementRadiusSize) >= stimulus.minRadius
                        radius1 = max(radius,stimulus.minRadius);
                        radius2 = min(radius+stimulus.elementRadiusSize,stimulus.maxRadius);
                        
                        %calculate in polar angle coordinates the corners of this quad
                        r = [radius1 radius1 radius2 radius2];
                        a = [angle angle+stimulus.elementAngleSize angle+stimulus.elementAngleSize angle];
                        
                        %convert into rectilinear coordinates and save in array
                        stimulus.x{phaseNum}(:,end+1) = r.*cos(d2r(a));
                        stimulus.y{phaseNum}(:,end+1) = r.*sin(d2r(a));
                        
                        % also calculate what contrast we want
                        stimulus.c{phaseNum}(:,end+1) = [1 1 1]*(isodd(radiusNum+isodd(angleNum)));
                        
                    end
                end
            end
            disppercent(length(allPhases1));
            disppercent(inf);
            stimulus.n = length(allPhases1);
            stimulus.phaseNum = 1;
            stimulus.flickerNum = 1;
        end
%     end
else
    disp(sprintf('(mglRetinotopy) Using precomputed stimulus pattern'));
end

%remember these parameters, so that we can know whether we
%need to recompute
stimulus.last.elementRadiusSize = stimulus.elementRadiusSize;
stimulus.last.elementAngleSize = stimulus.elementAngleSize;
stimulus.last.elementRadialVelocity = stimulus.elementRadialVelocity;
stimulus.last.maxRadius = stimulus.maxRadius;
stimulus.last.minRadius = stimulus.minRadius;
stimulus.last.elementVelocity = stimulus.elementVelocity;

% new we calculate the masks that cover the stimulus
% angles = (0:stimulus.wedgeStepSize:(360-stimulus.wedgeStepSize))+90+stimulus.wedgeAngle/2;
angles = task{2}{1}.randVars.myRandomloc;

% create masks for wedges (changes each trial)
for angleNum = 1:length(angles)
    angle = angles(angleNum);
    
    % init the wedge mask values
    stimulus.maskWedgeX{angleNum} = [];
    stimulus.maskWedgeY{angleNum} = [];
    
    % create a polygon that spares the wedge that we want
    % start in the center, compute it in radial coordinates
    r = 0;a = 0;
    
    % and go around the angles except for the wedge we want
    %vertexAngle = angle:(angle+360-stimulus.wedgeAngle);
    for vertexAngle = angle+stimulus.wedgeAngle/2:(angle+360-stimulus.wedgeAngle/2);
        r(end+1) = stimulus.maxRadius+1;
        a(end+1) = vertexAngle;
    end
    
    % and end up in the center
    r(end+1) = 0;
    a(end+1) = 0;
    
    % now convert to cartesian coordinates
    stimulus.maskWedgeX{angleNum}(:,end+1) = r.*cos(d2r(a));
    stimulus.maskWedgeY{angleNum}(:,end+1) = r.*sin(d2r(a));
end
stimulus.wedgeN = length(angles);
stimulus.currentMask = 0;
stimulus.currentFlicker = 0;

%update stimulus (each frame)
function [stimulus] = updateWedge(stimulus,myscreen)

%CONTRAST
%
%calculate max and min wedge element luminance for current contrast
%with min and max equidistant from the background luminance (grey (0.5)).
%e.g., if contrast is 10%, max and min luminance deviate each 5% from the
%grey background, thus 5+5 = 10% difference (contrast) between max
%and min.
LumMax = stimulus.dots.LumMax;
LumMin = stimulus.dots.LumMin;

%assign a luminance that produces the desired contrast
%to each element of the wedge
stimulus.cthisSeg = stimulus.c;

%update the phase of the sliding wedges (if1um,stimulus.n);
if isodd(stimulus.currentFlicker)
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 1) = LumMax;
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 0) = LumMin;
    %draw the whole stimulus pattern in black and white
    mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.cthisSeg{stimulus.phaseNum},0);
else
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 1) = LumMin;
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 0) = LumMax;    
    %draw the whole stimulus pattern in white and black (colors are reversed)
    mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.cthisSeg{stimulus.phaseNum},0);
end
%mask out to get a wedge (gray 0.5)
mglPolygon(stimulus.maskWedgeX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskWedgeY{stimulus.currentMask}+stimulus.yOffset,stimulus.background);





%convert from mouse's to aperture's coord
function mouseinfo = pos2ap(myscreen,fixStimulus,mouseinfo,stimulus)

%and the mouse coordinates relative to the center of the target screen
mouseinfo.x.pix = mouseinfo.x.pix-fixStimulus.pospix(1); %(in pix)
mouseinfo.y.pix = mouseinfo.y.pix-fixStimulus.pospix(2);

%convert pixels 2 angles (mglLines works with angles)
mouseinfo.x.angle.screen = mouseinfo.x.pix*mglGetParam('xPixelsToDevice');
mouseinfo.y.angle.screen = mouseinfo.y.pix*mglGetParam('yPixelsToDevice');

%calculate the transformation parameter that position the mouse on the aperture
%calculate length of the arrow (pythagoras theorem)
arrowsz2 = mouseinfo.x.angle.screen^2 + mouseinfo.y.angle.screen^2;%(angle)

%calculate transformation parameter
transpara = stimulus.dots.rmax^2/arrowsz2;

%transform actual coordinates to on-aperture coordinates.
mouseinfo.x.angle.onap=sqrt(mouseinfo.x.angle.screen^2*transpara)*sign(mouseinfo.x.angle.screen);
mouseinfo.y.angle.onap=sqrt(mouseinfo.y.angle.screen^2*transpara)*sign(mouseinfo.y.angle.screen);

%convert from PowerMate's to aperture's coord
function mouseinfo = PowerMate2ap(mouseinfo,stimulus,task)
global Inicoord;
%Set speed, initial position and size of arrow.
%wheelspeed=0.008;
wheelspeed=0.032;
%fixed
%initAngle=(225+270)/2*pi/180;
%random
initAngle.rad=task.thistrial.initAngledeg*pi/180;

r=stimulus.dots.rmax;

%Calculate arrow's coordinates on aperture
mouseinfo.x.angle.onap= r * cos(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
mouseinfo.y.angle.onap= r * sin(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap));

%convert from polar to cartesian coord
function coord = polar2cartesian(theta,r)
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

%set new Gamma table for better luminance resolution
function setGammaTableForMaxContrast(maxContrast)
%sets the gamma table so that we can have
%finest possible control over the stimulus luminance thus contrast.
global stimulus;

%set reserved colors RGB triplets
%at the bottom of the table (first rows)
gammaTable(1:size(stimulus.colors.resvdTriplet,1),1:size(stimulus.colors.resvdTriplet,2)) = stimulus.colors.resvdTriplet;

%set gamma table
if maxContrast > 0
  
  %fill up the rest of the gamma table
  cmax = 0.5 + maxContrast/2; 
  cmin = 0.5 - maxContrast/2;
  luminanceVals = cmin:((cmax - cmin)/(stimulus.colors.nPossCol-1)):cmax;

  %now get the linearized range
  %The gamma table has 256 possible luminance values from 0:1/255:1
  %We interpolate linearly the npossible  red/green/blue values at luminance values = luminanceVals
  redLinearized   = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
  greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
  blueLinearized  = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
  
  %add these values to the table
  gammaTable((stimulus.colors.minStimIx:stimulus.colors.mxSIx)+1,:) = [redLinearized;greenLinearized;blueLinearized]';
end

%set gamma table
mglSetGammaTable(gammaTable);

%remember what current max contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%get contrast index in new gammaTable
function [actualCon,LumMaxIx,LumMinIx] = getNormIdxContrast(curContrast,gammaTableNew)

global stimulus
% mglClose
% keyboard

%new table gray
gray = (gammaTableNew(stimulus.colors.midStimIx));

%max and min table luminances needed to implement current contrast
LumMax = gray + curContrast/2;
LumMin = gray - curContrast/2;

%indices in the table
[~,ixMax] = min(abs(bsxfun(@minus,gammaTableNew,[LumMax LumMax LumMax])));
[~,ixMin] = min(abs(bsxfun(@minus,gammaTableNew,[LumMin LumMin LumMin])));

%check unique indices
if length(unique(ixMax))==1 && length(unique(ixMin))==1
else
    mglClose
    keyboard
end

%actual contrast
actualCon = gammaTableNew(ixMax(1)) - gammaTableNew(ixMin(1));

%get index as normalized luminance positions
%ranging from 0 to 1. These are the luminances
%that we actually use in the code to display
%stimulus. When they are called, the hardware
%automatically check the associated gamma table 
%luminances
normIx = [0 [1:255]/255];
LumMaxIx = normIx(ixMax(1));
LumMinIx = normIx(ixMin(1));