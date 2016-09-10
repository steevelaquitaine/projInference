
%slMakeLocStimWdgAreaGrBg.m
%
%
%author: steeve laquitaine
%  date: 151129
% usage:
%
%       hi-pass
%       slMakeLocStimWdgAreaGrBg('displayName=full','sigmaw=30','loadPrecompStim');%better ! (more gradient ?)
%
%       lo-pass
%       slMakeLocStimWdgAreaGrBg('displayName=test','sigmaw=3')
%
%
%options
%
%       'loadPrecompStim': load precomputed stimulus from
%                          'stim_slMakeLocStimWdgAreaGrBg.mat' file
%
%Description:
%
%       location stimulus contrast is such that within wedge area
%       region lighter than the mean luminance within wedge area
%       see an increase in their luminance by contrast/2 and darker
%       one a decrease by contrast/2.
%       gray bkg.
%       can use 3 con ranging from 0 to 1 (0 to 100%)


function myscreen = slMakeLocStimWdgAreaGrBg(varargin)

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
task{2}{1}.segmin = [1 2];
task{2}{1}.segmax = [1 2];
params.task.parameter.loc.series = repmat([5:10:355]',3,1);
numtrials = numel(params.task.parameter.loc.series);
params.task.parameter.loc.trial = 1:numtrials;
params.task.parameter.loc.con   = repmat([1 0.12 0.1],36,1);%repmat([0.0009 0.0018 0.0036],36,1);
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

%screen size settings
stimulus.xOffset = xOffset;
stimulus.yOffset = yOffset;

%init stim
myscreen = initStimulus('stimulus',myscreen);

%------------------------- fixation --------------------------
myscreen.fixCr   = stimulus.colors.resvdNormIx(2);
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
% if ieNotDefined('params')
%     myTrialRand = 1 : numel(params.task.parameter.loc.series);
% end

myTrialRand = 1 : numel(params.task.parameter.loc.series);

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

% mglClose; keyboard
%case load precomputed stimulus
if slIsInput(varargin,'loadPrecompStim')
    load('stim_slMakeLocStimWdgAreaGrBg')
    slPrintfStr('slMakeLocStimWdgAreaGrBg','Loading precomp stimuli')    
    %input are indices between 5 and 255
    for i = 1 : numtrials
        task{2}{1}.randVars.tex{i} = mglCreateTexture(flipud(scaledbI255{i}));
    end
else    
    %Compute stim
    slPrintfStr('slMakeLocStimWdgAreaGrBg','Computing stimuli')
    for i = 1 : numtrials
        %status
        sldispPerc('slMakeLocStim',i,numtrials);
        %create stimulus here to add to noisy background texture
        loc = task{2}{1}.randVars.myRandomloc(i);
        tp = rand(ImgH,ImgW)>0.5; %rand pix img
        angMask = slGraphWedges2(ImgH,ImgW,400,loc,5);
        %fft wedge mask
        angMaskft = slFourierFilt(tp,sigmaw,'gausWeights'); %low pass filter (FFT)
        if isempty(angMaskft)
            slPrintfStr('slMakeLocStim','Input weights... For now weights can only be "gausWeights"')
            mglClose; keyboard
        end
        %rescale bkg luminances between 0 and 1
        %note: fourier filtering biases white and black
        %luminances toward gray luminances
        mx2 = max(angMaskft(:)); mn2 = min(angMaskft(:));
        scaledFiltMask = (angMaskft-mn2)./(mx2-mn2);
        bkg = ones(ImgH,ImgW)*stimulus.colors.grayColor; %gray
        bkg(angMask==1) = scaledFiltMask(angMask==1);  %cloudy mask
        %smooth it
        %bkgSmth = slFourierFilt(bkg,10,'gausWeights'); %low pass filter (FFT)
        bkgSmth =bkg;
        %zscore to get equal proba of
        %of darker and lighter luminance relative
        %to gray background
        bkgSmthZ = zeros(ImgH,ImgW);
        bkgSmthZ = (bkgSmth - mean(bkgSmth(:)))/std(bkgSmth(:)); %0 mean 1std
        %scale between 0 and 255 with gray background mean
        b = bkgSmthZ + stimulus.colors.midStimIx;
        %scale lighter areas
        %betw. 0.5 and 1
        li=b>mean(b(:));%stimulus.colors.midStimIx;
        mxli = max(b(li));
        %scaledbLi = 0.5*(b(li)-stimulus.colors.midStimIx)./(mxli - stimulus.colors.midStimIx);
        scaledbLi = 0.5*(b(li) - mean(b(:)))./(mxli - mean(b(:)));
        %scale darker ones
        %between 0 and 0.5
        dk=b<mean(b(:));%stimulus.colors.midStimIx;
        mndk = min(b(dk));
        %scaledbDk = 0.5*(b(dk)-stimulus.colors.midStimIx)./(mndk - stimulus.colors.midStimIx);
        scaledbDk = 0.5*(b(dk) - mean(b(:)))./(mndk - mean(b(:)));
        scaledAll = ones(ImgH,ImgW)*0.5; %mean is 0.5
        scaledAll(li) = scaledbLi + 0.5;
        scaledAll(dk) = 0.5 - scaledbDk;
        scaledbI255{i} = scaledAll*(stimulus.colors.nPossCol-1)+stimulus.colors.nresvd+1;
        
        %input are indices between 5 and 255
        task{2}{1}.randVars.tex{i} = mglCreateTexture(flipud(scaledbI255{i}));
    end
end
save('stim_slMakeLocStimWdgAreaGrBg','scaledbI255')

% mglClose; keyboard
% task{2}{1}.randVars.actualCon = actualCon;
slPrintfStr('slMakeLocStim',' done.')


%--------------------------- start ------------------------------------
%Init task
[task{2}{1},myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);

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
    %     stimulus.loc = task.thistrial.myRandomloc;
end
stimulus.loc = task.thistrial.myRandomloc;

%------------------------
%run at each screen frame
%------------------------
function [task,myscreen] = updateScreenCallback(task, myscreen)
global stimulus;
mglClearScreen(stimulus.colors.grayColor)
%need to be re-set after eyecalib
%setGammaTableForMaxContrast(stimulus.maxCon);
setGammaTableForMaxContrast(task.thistrial.myRandomCon);

%location
if (task.thistrial.thisseg == 2)
    %coord centered at screen center
    mglBltTexture(task.thistrial.tex,[0 0]);
end

%------------------------
%get response
%------------------------
function [task,myscreen] = getResponseCallBack(task,myscreen)

%When subject presses key, stimulus again
task = jumpSegment(task,1);



%set new Gamma table for better luminance resolution
function setGammaTableForMaxContrast(maxContrast)
%sets the gamma table so that we can have
%finest possible control over the stimulus luminance thus contrast.
global stimulus;

%set reserved colors RGB triplets
%at the bottom of the table (first rows)
gammaTable(1:size(stimulus.colors.resvdTriplet,1),1:size(stimulus.colors.resvdTriplet,2)) = stimulus.colors.resvdTriplet;
% mglClose; keyboard

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
    gammaTable((stimulus.colors.minStimIx:stimulus.colors.mxSIx),:) = [redLinearized;greenLinearized;blueLinearized]';
end

%set gamma table
mglSetGammaTable(gammaTable);

%remember what current max contrast is that we can display
stimulus.currentMaxContrast = maxContrast;





