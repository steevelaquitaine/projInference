
%slMakeLocStimWdgArea.m
%
%
%author: steeve laquitaine
%  date: 151112
% usage:
%       
%       hi-pass        
%       slMakeLocStimWdgArea('displayName=test','sigmaw=30');better ! (more gradient ?)
%
%       lo-pass
%       slMakeLocStimWdgArea('displayName=test','sigmaw=2')
%
%Description:
%           
%       My favorite so far: 
%       location stimulus contrast is such that within wedge area
%       region lighter than the mean luminance within wedge area
%       see an increase in their luminance by contrast/2 and darker
%       one a decrease by contrast/2.


function myscreen = slMakeLocStimWdgArea(varargin)

%get args
eval(evalargs(varargin,0,0,{'displayName','sigmaw'}));
myscreen.displayname = displayName;
myscreen = initScreen(myscreen);

%check args
if ieNotDefined('sigmaw')
    slPrintfStr('slMakeLocStim',' ! WARNING ! Please define sigmaw (e.g., "sigmaw=3") ...')
    mglClose
    keyboard
end

%-----------------------------init Task----------------------------
task{2}{1}.segmin = [1 .3];
task{2}{1}.segmax = [1 .3];
params.task.parameter.loc.series = repmat([5:10:355]',3,1);
numtrials = numel(params.task.parameter.loc.series);
params.task.parameter.loc.trial = 1:numtrials;
params.task.parameter.loc.con   = repmat([0.1 0.02 0.01],36,1);%repmat([0.0009 0.0018 0.0036],36,1);
params.task.parameter.loc.mean  = 225;
params.task.parameter.loc.modes = 225;
params.task.parameter.loc.std   = inf;
params.task.parameter.loc.strength = 0;
params.task.parameter.loc.sample.degree = unique(params.task.parameter.loc.series);


%-------------------------- Stimulus -------------------------------
global stimulus;

%init gamma table
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%Colors & contrasts
stimulus.colors.resvdTriplet = [0 0 0; 1 1 1; 0 1 0; 1 0 0];
stimulus.colors.nresvd = size(stimulus.colors.resvdTriplet,1);
stimulus.colors.resvdIx = [1:4];
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
if ieNotDefined('eltSize'),eltSize = [];end
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
stimulus.minRad = 0;
stimulus.maxRad = stimulus.dots.rmax;
stimulus.sizeStim = sizeStim;

%Wedge elements
stimulus.eltAngSz      = 5; %size in degrees
stimulus.eltRadSize    = 2; %radial length
stimulus.eltRadialVelocity = 0; %radial speed (0 for static)
stimulus.eltVelocity   = 6; %elt velocity
stimulus.FlickerFreq   = 8;%8 %flickering frequency in hertz

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
myscreen.fixRad     = 0.2;
myscreen.fixPosX    = 0;
myscreen.fixPosY    = 0;
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


%-------------------------- Noisy background -------------------------------
%get dim (pix)
ImgW = myscreen.screenWidth;
ImgH = myscreen.screenHeight;

%precompute a background
%for each trial
slPrintfStr('slMakeLocStim',' Creating noisy Backg...')
stimulus.maxCon = max(task{2}{1}.randVars.myRandomCon);

%test----------
numtrials = 10;
%--------------

for i = 1 : numtrials    
    %status
    sldispPerc('slMakeLocStim',i,numtrials);    
    %make background
    bkg = rand(ImgH,ImgW)>0.5; %rand pix img
    bkgf = slFourierFilt(bkg,sigmaw,'gausWeights'); %low pass filter (FFT)
    if isempty(bkgf)
        slPrintfStr('slMakeLocStim','Input weights... For now weights can only be "gausWeights"')
        mglClose; keyboard
    end
    %rescale bkg luminances between 0 and 1
    %note: fourier filtering biases white and black
    %luminances toward gray luminances
    mx = max(bkgf(:)); mn = min(bkgf(:));
    bkgfsc = (bkgf-mn)./(mx-mn);        
    %scale between min and max luminance indices
    %between nres + 1 and 255
    bkg2 = bkgfsc*(stimulus.colors.nPossCol-1)+stimulus.colors.nresvd+1;        
    %create stimulus here to add to noisy background texture
    loc = params.task.parameter.loc.series(i);
    angMask = slGraphWedges2(ImgH,ImgW,200,loc,10);    
    %enhance contrast of a wedge region
    [actualCon(i),angMaskIx255] = getNormIdxContrast(task{2}{1}.randVars.myRandomCon(i),gammaTableNew,bkg2,angMask);
    bkg2(angMask==1) = angMaskIx255;    
    %input are indices between 5 and 255
    task{2}{1}.randVars.tex{i} = mglCreateTexture(bkg2);    
end
task{2}{1}.randVars.actualCon = actualCon;
slPrintfStr('slMakeLocStim',' done.')


%--------------------------- start ------------------------------------
%Init task
[task{2}{1},myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);

%Init stim
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

%------------------------
%run at each trial segment
%------------------------
function [task,myscreen] = startSegmentCallback(task,myscreen)
global stimulus;
global fixStimulus;
global Inicoord;

%Stimulus
if (task.thistrial.thisseg == 2)
    
    %------ wedge stim ---------
    %update location strength each trial
    stimulus.dots.contrast = task.thistrial.actualCon;
%     stimulus.dots.LumMax = task.thistrial.LumMax;
%     stimulus.dots.LumMin = task.thistrial.LumMin;
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

%------------------------
%run at each screen frame
%------------------------
function [task,myscreen] = updateScreenCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.colors.grayColor)
%need to be re-set after eyecalib
setGammaTableForMaxContrast(stimulus.maxCon);

%location
if (task.thistrial.thisseg == 2)    
    
    %------ Noisy bkg + stim ---------
    %draw bkg
    %coord centered at screen center
    mglBltTexture(task.thistrial.tex,[0 0]);

    %update frame count each frame
    stimulus.countFrame = stimulus.countFrame + 1;
    %update stimulus only when number of frames satisfying our flickering
    %frequency is reached otherwise stimulus remains frozen
    if SLisMultiple(stimulus.countFrame,stimulus.nFrameStill);
        stimulus.currentFlicker = stimulus.currentFlicker + 1;
    end
end

%------------------------
%get response
%------------------------
function [task,myscreen] = getResponseCallBack(task,myscreen)

%When subject presses key, stimulus again
task = jumpSegment(task,1);



%init stim
function stimulus = initWedge(stimulus,myscreen,task)

%round stimulus radius to nearest quarter of a degree,
%this reduces some edge effects
stimulus.maxRad = floor(stimulus.maxRad/.25)*.25;
disp(sprintf('Stimulus radius = [%0.2f %0.2f] degrees',stimulus.minRad,stimulus.maxRad));

%calculate params
%size of wedges
stimulus.wedgeAngle = 360*stimulus.sizeStim;

%how much to step the wedge angle by
%stimulus.wedgeStepSize = 360/stimulus.stepsPerCycle;
sizeStim = stimulus.sizeStim;

%we only need to recompute the mglQuad points of the elts if something has
%changed in the stimulus. This is for the radial elt pattern
if ~isfield(stimulus,'last') || ~isfield(stimulus,'x') || ...
        (stimulus.eltAngSz ~= stimulus.last.eltAngSz) || ...
        (stimulus.eltRadSize ~= stimulus.last.eltRadSize) || ...
        (stimulus.maxRad ~= stimulus.last.maxRad) || ...
        (stimulus.minRad ~= stimulus.last.minRad)
    
    %all the angles that the elts will be made up of
    allAngles = (0:stimulus.eltAngSz:(360-stimulus.eltAngSz));
    
    %all the phases. The phase refers to the radial position of the
    %black and white pattern.
    %static stimulus
    allPhases1 = 1;
    
    %# of flickering
    allFlickers = task{2}{1}.segmax*myscreen.framesPerSecond;
    disppercent(-inf,'(slMakeLocStim) Calculating coordinates of elements in stimulus pattern');
    phaseNum = 1;
    stimulus.x{phaseNum} = [];
    stimulus.y{phaseNum} = [];
    stimulus.c{phaseNum} = [];
    
    %for each angle that the elt contains
    for angleNum = 1:length(allAngles)
        
        %get the angle
        angle = allAngles(angleNum);
        
        %choose which phase we are going to be
        %thisMinRad = stimulus.minRad - allPhases1(phaseNum);
        thisMinRad = stimulus.minRad - allPhases1;
        
        %all the radiuses
        allRad = thisMinRad : stimulus.eltRadSize : stimulus.maxRad;
                
        %now create all the quads for this wedge
        for radiusNum = 1 : length(allRad)
            radius = allRad(radiusNum);
            
            if (radius+stimulus.eltRadSize) >= stimulus.minRad
                radius1 = max(radius,stimulus.minRad);
                radius2 = min(radius+stimulus.eltRadSize,stimulus.maxRad);
                
                %calculate in polar angle coordinates the corners of this quad
                r = [radius1 radius1 radius2 radius2];
                a = [angle angle+stimulus.eltAngSz angle+stimulus.eltAngSz angle];
                
                %convert into rectilinear coordinates and save in array
                stimulus.x{phaseNum}(:,end+1) = r.*cos(d2r(a));
                stimulus.y{phaseNum}(:,end+1) = r.*sin(d2r(a));
                
                %also calculate what contrast we want
                stimulus.c{phaseNum}(:,end+1) = [1 1 1]*(isodd(radiusNum+isodd(angleNum)));
            end
        end
    end
    disppercent(length(allPhases1));
    disppercent(inf);
    stimulus.n = length(allPhases1);
    stimulus.phaseNum = 1;
    stimulus.flickerNum = 1;
else
    disp(sprintf('(mglRetinotopy) Using precomputed stimulus pattern'));
end

%remember these parameters, so that we can know whether we
%need to recompute
stimulus.last.eltRadSize = stimulus.eltRadSize;
stimulus.last.eltAngSz = stimulus.eltAngSz;
stimulus.last.eltRadialVelocity = stimulus.eltRadialVelocity;
stimulus.last.maxRad = stimulus.maxRad;
stimulus.last.minRad = stimulus.minRad;
stimulus.last.eltVelocity = stimulus.eltVelocity;

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
        r(end+1) = stimulus.maxRad+1;
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
%calculate max and min wedge elt luminance for current contrast
%with min and max equidistant from the background luminance (grey (0.5)).
%e.g., if contrast is 10%, max and min luminance deviate each 5% from the
%grey background, thus 5+5 = 10% difference (contrast) between max
%and min.
LumMax = stimulus.dots.LumMax;
LumMin = stimulus.dots.LumMin;

%assign a luminance that produces the desired contrast
%to each elt of the wedge
stimulus.cthisSeg = stimulus.c;

%update phase of sliding wedges
if isodd(stimulus.currentFlicker)
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 1) = LumMax;
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 0) = LumMin;
    %draw whole pattern in black and white
    mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.cthisSeg{stimulus.phaseNum},0);
else
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 1) = LumMin;
    stimulus.cthisSeg{stimulus.phaseNum}(stimulus.c{stimulus.phaseNum} == 0) = LumMax;
    %draw the whole stimulus pattern in white and black (colors are reversed)
    mglQuad(stimulus.x{stimulus.phaseNum}+stimulus.xOffset,stimulus.y{stimulus.phaseNum}+stimulus.yOffset,stimulus.cthisSeg{stimulus.phaseNum},0);
end

%mask out to get a wedge (gray 0.5)
mglPolygon(stimulus.maskWedgeX{stimulus.currentMask}+stimulus.xOffset,stimulus.maskWedgeY{stimulus.currentMask}+stimulus.yOffset,stimulus.background);
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
%for each element of the location mask
function [actualCon,angMaskIx255] = getNormIdxContrast(curContrast,gammaTableNew,backLum,angMask)

global stimulus

%new table gray
gray = (gammaTableNew(stimulus.colors.midStimIx));

%max and min table luminances needed to implement current contrast
angMaskIx255 = round(backLum(angMask==1));
angMaskIx1 = gammaTableNew(angMaskIx255,1);
% LumMax = angMaskIx1+0.4; 
% LumMin = angMaskIx1-0.4; 
LumMax = angMaskIx1 + curContrast/2;
LumMin = angMaskIx1 - curContrast/2;

%possible luminance for contrast
%in the gamma table
gammaPoss = gammaTableNew(stimulus.colors.minStimIx+1:end,:);

%indices in the table
for i=1:length(LumMax)
    [~,ixMax(i,:)] = min(abs(bsxfun(@minus,gammaPoss,[LumMax(i) LumMax(i) LumMax(i)])));
    [~,ixMin(i,:)] = min(abs(bsxfun(@minus,gammaPoss,[LumMin(i) LumMin(i) LumMin(i)])));
end

%express ix between 5 and 255
%for mglCreateTexture
ixMax = ixMax(:,1) + stimulus.colors.nresvd;
ixMin = ixMin(:,1) + stimulus.colors.nresvd;

%increase contrast between lighter and darker regions
%than gray within angle mask (the location region)
midIxAngMask = round(mean(angMaskIx255));%(max(ixMax)-min(ixMin))/2;
lighterPixs = angMaskIx255>midIxAngMask;%stimulus.colors.midStimIx;
darkerPixs = angMaskIx255<midIxAngMask;%stimulus.colors.midStimIx;
angMaskIx255(lighterPixs) = ixMax(lighterPixs);
angMaskIx255(darkerPixs) = ixMin(darkerPixs);

%actual average contrast within region
AvLumMax = mean(gammaTableNew(angMaskIx255(lighterPixs)));
AvLumMin = mean(gammaTableNew(angMaskIx255(darkerPixs)));
actualCon = AvLumMax - AvLumMin;








