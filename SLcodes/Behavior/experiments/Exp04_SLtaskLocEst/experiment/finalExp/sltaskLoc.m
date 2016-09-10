
%sltaskLoc.m
%
%
% author: steeve laquitaine
%   date: 151202
%purpose: noisy location stimulus task
% status: complete
%   note: stimuli variable can reach 3GB for large screen!
%         making stim takes 1 min
%
%  usage:
%
%       sltaskLoc('displayName=Test','loadStim','p=10');
%       sltaskLoc('displayName=Test','computeStim','p=10');
%       sltaskLoc('displayName=VPixx','computeStim','p=80');
%       sltaskLoc('displayName=VPixx','loadStim','p=80');
%
%
%
%Mandatory inputs:
%           'displayName': name of display screen
%                          e.g., 'displayName=steeveScreen'
%                     'p': 80 (40/20 or 10) load params file created by
%                     SLinitRunExpUniPriorLoc.m
%                          e.g., 'p=80'
%       'loadPrecompStim': load precomputed stimulus from
%                          'stim_slMakeLocStimWdgAreaGrBg.mat' file
%
%
%        Option: 'responseDevice':'responseDevice=1'(0) for powermate
%        (mouse) 'responseObject':'responseObject=arrow'('dot')
%           'computeStim': make new stimuli (instead of 'loadPrecompStim')
%
%Description:
%
%       location stimulus contrast is such that within wedge area region
%       lighter than the mean luminance within wedge area see an increase
%       in their luminance by contrast/2 and darker one a decrease by
%       contrast/2. gray bkg. can use 3 con ranging from 0 to 1 (0 to 100%)


function myscreen = sltaskLoc(varargin)

%get inputs Default response device is powermate
displayName=[];p=[];responseDevice=[];responseObject=[];
getArgs(varargin,{'displayName','p=[]','responseDevice=1','responseObject=arrow'});

%get precomputed params from SLcodes library
if p==80
    params = 'steeve_exp12_metho_Pstd080_mean225_coh010012100_dir36_t100_075_034perCoh_131224.mat';
elseif p==40
    params = 'steeve_exp12_metho_Pstd040_mean225_coh010012100_dir36_t100_075_032perCoh_131224.mat';
elseif p==20
    params = 'steeve_exp12_metho_Pstd020_mean225_coh010012100_dir36_t103_074_029perCoh_131224.mat';
elseif p==10
    params = 'steeve_exp12_metho_Pstd010_mean225_con010012100_dir36_t107_074_033perCoh_131224.mat';
else
    fprintf('%s \n','(sltaskLoc) need to set p=80 or p=40 ...')
end


%---------- init screen ---------- displayName
myscreen.displayname = displayName;
%response device
myscreen.responseDevice = responseDevice;
%response screen
myscreen.responseObject= responseObject;
%init
myscreen = initScreen(myscreen);
%hide cursor
mglDisplayCursor(0);

%check screen is not flipped it changes locations
if ~sum(myscreen.flipHV)==0
    fprintf(['!! WARNING !! The screen is flipped.',...
        'This might flip the stimulus displayed location. \n'])
    fprintf('If you must flip the screen, check that the displayed location on screen matches the output location in the terminal. \n')
    mglClose; keyboard
end

%-------- init task params ------- 
%trial segments fix - stim - resp - conf
%- fb
task{2}{1}.segmin = [1 .02 5 .1 .1];
task{2}{1}.segmax = [1 .02 5 .1 .1];

%params
pf = load(params);
pf = pf.task.parameter.loc;

%randomize trials
nTrials = length(pf.series);
rndT = randperm(nTrials);

%params
task{2}{1}.randVars.myRandomloc = pf.series(rndT);
task{2}{1}.randVars.myRandomCon = pf.con(rndT);
task{2}{1}.randVars.myRandomCon(task{2}{1}.randVars.myRandomCon==0.12) = 0.156;

%debugging----------
%task{2}{1}.randVars.myRandomCon = ones(nTrials,1)*0.10;
%debugging----------
task{2}{1}.randVars.myMean      = pf.mean;
task{2}{1}.randVars.myModes     = 225;
task{2}{1}.randVars.myStrength  = pf.std;
task{2}{1}.randVars.p           = p;
task{2}{1}.randVars.myLocations = pf.sample;

%randomize mouse location at the start
nTrials = numel(task{2}{1}.randVars.myRandomloc);
task{2}{1}.randVars.initAngledeg = randi([0,359],[nTrials,1]);

%get response and reaction time response segment is set at "3" or "4" for
%mouse event and 1 for keyboard press. We use keyboard press "1" to enter
%response
task{2}{1}.getResponse = [0 0 1 0 0];

%init stored data location reported by subject in deg and in cartesian
task{2}{1}.randVars.calculated.prodcoor = [nan nan];
task{2}{1}.randVars.calculated.proddeg = nan;


%-------- init Stimulus --------
global stimulus;

%init gamma table
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%reserve 4 colors
stimulus.colors.resvdTriplet = [0 0 0; 1 1 1; 0 1 0; 1 0 0];
stimulus.colors.nresvd = size(stimulus.colors.resvdTriplet,1);

%# of stim. luminances (255 - # reserved)
maxIx = 255;
stimulus.colors.nPossCol = maxIx+1 - stimulus.colors.nresvd;

%min, mid and max index row stim colors in gammaTable (index values start
%from 0 based). e.g., min/mid/max/StimIx = 253/254/255;
stimulus.colors.minStimIx = stimulus.colors.nresvd+1;
stimulus.colors.midStimIx = stimulus.colors.minStimIx + floor(stimulus.colors.nPossCol/2);
stimulus.colors.mxSIx     = maxIx;
stimulus.colors.grayColor = stimulus.colors.midStimIx/maxIx;

%set the reserved colors as normalized luminance indices (always [0
%[1:255]/255])) Those indices are associated with each of the 255 rows in
%the gammaTable and with each of the corresponding luminances triplets for
%RGB. e.g., 0 is always associated with the gammaTable' first RGB luminance
%triplet (1st row) and 1 with the last (255th row)%
% e.g., the indices
%of 5 reserved colors positioned at the first rows of the gamma table will
%be [0:4]/255 Those are the indices called to display their associated
%luminances.
i = [0 [1:maxIx]/maxIx];
stimulus.colors.resvdNormIx = i(1:stimulus.colors.nresvd);

%convert distance in visual angle to pixels wedge radius in mm
stimulus.radiusdeg = 2.5; %deg
d2screen = myscreen.displayDistance*10; %mm
stimulus.radiusmm = slGetVisA2mm(d2screen,stimulus.radiusdeg); %mm
PixBymm = slGetScreenPixbyMM(myscreen.screenNumber);
stimulus.radiuspix = stimulus.radiusmm*PixBymm; %pixels

%precompute n trials stimuli based on params
if slIsInput(varargin,'computeStim')
    slPrintfStr('sltaskLoc',' Computing stimuli. Takes time...')
    t1=tic;
    mytex = computeStim(task,myscreen,stimulus,nTrials);
    task{2}{1}.randVars.tex = mytex;
    tend=toc(t1);
    slPrintfStr('sltaskLoc',['Took ' num2str(tend) 'sec'] )
elseif slIsInput(varargin,'loadStim')
    %check stim file exists
    ise = dir(['stim_sltaskLoc' num2str(p) '.mat']);
    if isempty(ise);
        fprintf('%s \n',['(sltaskLoc) You can not call ',...
            '"loadStim" because "stim_sltaskLoc" was not ',...
            'found in cd. You can either get the file ',...
            'somewhere or call "computeStim" to ',...
            'recalculate the stim'])
        mglClose; keyboard
    end
    %load stim
    load(['stim_sltaskLoc' num2str(p)]);
    %re-create texture cannot put in parfor because the code can't find the
    %opened screen
    for i = 1 : nTrials
        %status
        sldispPerc('slMakeLocStim',i,nTrials);
        %image must be flipped (reversed) in the height dimension to be
        %displayed properly because matlab matrices increase y indices from
        %top to bottom but openGL read in reverse direction (?)
        mytex{i} = mglCreateTexture(flip(scaledbI255{i},3));
    end
    task{2}{1}.randVars.tex = mytex;
end
nTrials = length(task{2}{1}.randVars.tex);
%reset # of trials to # of stimulus
task{2}{1}.numTrials = nTrials;

%init stim
myscreen = initStimulus('stimulus',myscreen);

%-------- init fixation -------
myscreen.fixCr   = stimulus.colors.resvdNormIx(1);
myscreen.fixRad  = 0.2;
myscreen.fixPosX = 0;
myscreen.fixPosY = 0;
[task{1},myscreen] = fixationLoc(myscreen);

%--------------------------- start ------------------------------------
%Init task
[task{2}{1},myscreen] = initTask(task{2}{1},myscreen,@startSegmentCallback,@updateScreenCallback,@getResponseCallBack);

%run eye calibration
myscreen = eyeCalibDisp(myscreen);

%init output in terminal
fprintf('%s \n','---------------- Task status ----------')
fprintf('%s   %s    %s  %s  %s  \n','Trial','Con','Loc','Resp','RT')

%Main display loop
phaseNum = 1;
while (phaseNum <= length(task{2})) && ~myscreen.userHitEsc
    
    %update task
    [task{2},myscreen,phaseNum] = updateTask(task{2},myscreen,phaseNum);
    
    %update the fixation task
    [task{1},myscreen] = updateTask(task{1},myscreen,1);
    myscreen = tickScreen(myscreen,task);
end

%End
myscreen = endTask(myscreen,task);
clear global stimulus;
clear global fixStimulus;



%------------------ run at each trial segment -------------------------
function [task,myscreen] = startSegmentCallback(task,myscreen)
global stimulus;
global fixStimulus
global Inicoord;

stimulus.loc = task.thistrial.myRandomloc;

%response
if (task.thistrial.thisseg == 3)
    %Position mouse randomly Calculate angle coordinates on aperture
    %relative to center
    coord = polar2cartesian(task.thistrial.initAngledeg,stimulus.radiusdeg);
    Inicoord.x.angle.onap = coord.x;
    Inicoord.y.angle.onap = coord.y;
    %Calculate pixels coordinates on aperture relative to center
    Inicoord.x.pix.onap = Inicoord.x.angle.onap * mglGetParam('xDeviceToPixels');
    Inicoord.y.pix.onap = Inicoord.y.angle.onap * mglGetParam('yDeviceToPixels');
    %Calculate coordinates on aperture relative to screen's root.
    Inicoord2root.x.pix.onap = Inicoord.x.pix.onap + fixStimulus.pospix(1);
    Inicoord2root.y.pix.onap = Inicoord.y.pix.onap + fixStimulus.pospix(2);
    %Position mouse
    mglSetMousePosition(Inicoord2root.x.pix.onap,Inicoord2root.y.pix.onap, myscreen.screenNumber);
end

%confirmation
if (task.thistrial.thisseg == 4)
    %Get mouse position
    mi = mglGetMouse(myscreen.screenNumber); %pix positions
    mouseinfo.x.pix = mi.x;
    mouseinfo.y.pix = mi.y;
    %Check if subject confirmed his choice: ("space bar is down") if you
    %want to you use mouse click to enter choice mouseinfo.buttons =
    %mi.buttons; if you want to you use keyboard button '1' to enter choice
    mouseinfo.buttons = mglGetKeys(19);
    %place mouse on aperture
    mouseinfo = pos2ap(fixStimulus,mouseinfo,stimulus);
    %interface PowerMate
    if myscreen.responseDevice == 1
        mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
    end
    %get response (cart)
    task.thistrial.prodcoor = [mouseinfo.x.angle.onap mouseinfo.y.angle.onap];
    %get response (deg)
    [~,task.thistrial.proddeg] = SLcart2polar(task.thistrial.prodcoor);
    %get reaction time
    if ~isfield(task.thistrial,'reactionTime')
        task.thistrial.reactionTime = NaN;
    end
    %print trial # with key and direction choosen.
    fprintf(' %i    %1.03f   %i   %i  %.1f \n',task.trialnum,...
        task.thistrial.myRandomCon,...
        task.thistrial.myRandomloc,...
        round(task.thistrial.proddeg),...
        task.thistrial.reactionTime)
    %then jump to feedback
    task = jumpSegment(task,4);
end

%------------------ run at each screen frame -------------------
function [task,myscreen] = updateScreenCallback(task,myscreen)
global stimulus;
global fixStimulus;

%Background luminance is always the mid index of the gamma table (gray)
mglClearScreen(stimulus.colors.grayColor)

%need to be re-set after eyecalib this sets the max/min stimulus luminances
%at each trial (thus contrast).
setGammaTableForMaxContrast(task.thistrial.myRandomCon);
%setGammaTableForMaxContrast(1);

%display stim. Stimulus max|min luminances is always this trial max|min
%gamma table luminance. The table max|min luminances change at each trial
%to implement desired contrast
if task.thistrial.thisseg == 2
    mglBltTexture(task.thistrial.tex,[0 0]);
end

%Response
if (task.thistrial.thisseg == 3)
    %place mouse on aperture
    mi = mglGetMouse(myscreen.screenNumber); %(pixels)
    mouseinfo.x.pix = mi.x;
    mouseinfo.y.pix = mi.y;
    mouseinfo = pos2ap(fixStimulus,mouseinfo,stimulus);
    %interface PowerMate
    if myscreen.responseDevice == 1
        mouseinfo = PowerMate2ap(mouseinfo,stimulus,task);
    end
    %draw a black arrow (radius)
    if strcmp(myscreen.responseObject,'arrow');
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            mouseinfo.x.angle.onap,...
            mouseinfo.y.angle.onap,2,stimulus.colors.resvdNormIx(1),1);
        %draw a white dot at the periphery (radius)
    elseif strcmp(myscreen.responseObject,'dot');
        mglPoints2(mouseinfo.x.angle.onap,mouseinfo.y.angle.onap,6,stimulus.colors.resvdNormIx(1))
    end
end

%confirmation
if (task.thistrial.thisseg == 4)
    %Position mouse on aperture deliver pixel positions, not angle
    mi = mglGetMouse(myscreen.screenNumber);
    mouseinfo.x.pix = mi.x;
    mouseinfo.y.pix = mi.y;
    mouseinfo = pos2ap(fixStimulus,mouseinfo,stimulus);
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
            mouseinfo.y.angle.onap,6,stimulus.colors.resvdNormIx(4))
    end
end

%feedback
if task.thistrial.thisseg == 5
    coord = polar2cartesian(task.thistrial.myRandomloc,stimulus.radiusdeg);
    feedbackcoord.x.angle.onap = coord.x;
    feedbackcoord.y.angle.onap = coord.y;
    %draw a green arrow (radius)
    if strcmp(myscreen.responseObject,'arrow')
        mglLines2(fixStimulus.pos(1),fixStimulus.pos(2),...
            feedbackcoord.x.angle.onap,...
            feedbackcoord.y.angle.onap,3,stimulus.colors.resvdNormIx(3),1);
        %draw a green dot at the periphery (radius)
    elseif strcmp(myscreen.responseObject,'dot');
        mglPoints2(feedbackcoord.x.angle.onap,feedbackcoord.y.angle.onap,8,stimulus.colors.resvdNormIx(3))
    end
end

%------------ get response ------------
function [task,myscreen] = getResponseCallBack(task,myscreen)

%When subject chose, confirm (jump to segment 4)
task=jumpSegment(task,4);




%make stimulus
function mytex = computeStim(task,myscreen,stimulus,nTrials)

%The resulting stimulus matrices always range between 5 (min reserved
%index) in the gamma table and 255 (max). The background is always the same
%mid index gray

%get task params
locAll = task{2}{1}.randVars.myRandomloc;
%get dim (pix)
ImgW = myscreen.screenWidth;
ImgH = myscreen.screenHeight;
%stim radius
radiuspix = stimulus.radiuspix;
%bkg color
grayColor = stimulus.colors.grayColor;
%mid stim index
midStimIx = stimulus.colors.midStimIx;
nPossCol = stimulus.colors.nPossCol;
nresvd = stimulus.colors.nresvd;
%create stim matrices
parfor i = 1 : nTrials
    %status
    sldispPerc('slMakeLocStim',i,nTrials);
    
    %get location
    loc = locAll(i);
    tp = rand(ImgH,ImgW)>0.5; %rand pix img
    angMask = slGraphWedges2(ImgH,ImgW,radiuspix,loc,5);
    
    %     %fft wedge mask low pass filter
    %     angMaskft = slFourierFilt(tp,30,'gausWeights');
    %
    %     %rescale bkg luminances between 0 and 1 because fourier filtering
    %     %biases white and black luminances toward gray luminances
    %     mx2 = max(angMaskft(:)); mn2 = min(angMaskft(:));
    %     scaledFiltMask = (angMaskft-mn2)./(mx2-mn2);
    %     bkg = ones(ImgH,ImgW)*grayColor;%gray
    %     bkg(angMask==1) = scaledFiltMask(angMask==1);   %cloudy mask
    
    %make sure scaledFiltMask range between 0 and 1
    scaledFiltMask = NaN;
    while max(scaledFiltMask(:))~=1 && min(scaledFiltMask(:))~=0
        %fft wedge mask low pass filter
        angMaskft = slFourierFilt(tp,30,'gausWeights');
        
        %rescale bkg luminances between 0 and 1 because fourier filtering
        %biases white and black luminances toward gray luminances
        mx2 = max(angMaskft(:)); mn2 = min(angMaskft(:));
        scaledFiltMask = (angMaskft-mn2)./(mx2-mn2);
    end      
    bkg = ones(ImgH,ImgW)*grayColor;%gray
    bkg(angMask==1) = scaledFiltMask(angMask==1);%cloudy mask
           
    %z-score to get equal proba of of darker and lighter luminance relative
    %to gray background
    bkgZ = zeros(ImgH,ImgW);
    bkgZ = (bkg - mean(bkg(:)))/std(bkg(:));
    
    %-----scale between 0 and 255 with gray mean bkg-----
    b = bkgZ + midStimIx;
    
    %scale lighter areas betw. 0.5 and 1
    li=b>mean(b(:));
    mxli = max(b(li));
    scaledbLi = 0.5*(b(li) - mean(b(:)))./(mxli - mean(b(:)));
    
    %scale darker ones
    dk=b<mean(b(:));
    mndk = min(b(dk));
    scaledbDk = 0.5*(b(dk) - mean(b(:)))./(mndk - mean(b(:)));
    scaledAll = ones(ImgH,ImgW)*0.5;%mean is 0.5
    scaledAll(li) = scaledbLi + 0.5;
    scaledAll(dk) = 0.5 - scaledbDk;
    
    %input are indices between 5 and 255.The file can reach 3GB so we uint8
    %to reduce the size format in RGBA x imagewidth x imageHeight necessary
    %to run with mglCreateTexture when using uint8 scaledAll must be
    %converted to width x height matrix by transpose
    scaledbI255WH = uint8((scaledAll')*(nPossCol-1)+nresvd+1);
    
    %RGBA x imagewidth x imageHeight
    scaledbI255tmRGB = repmat(scaledbI255WH,1,1,3);
    scaledbI255all = permute(scaledbI255tmRGB,[3 1 2]);
    scaledbI255all(4,:,:) = 255; %no alpha blending
    scaledbI255{i} = scaledbI255all;   
    
end

%create texture cannot put in parfor because the code can't find the opened
%screen
for i = 1 : length(scaledbI255)
    %status
    sldispPerc('slMakeLocStim',i,nTrials);
    %image must be flipped (reversed) in the height dimension to be
    %displayed properly because matlab matrices increase y indices from top
    %to bottom but openGL read in reverse direction (?)
    mytex{i} = mglCreateTexture(flip(scaledbI255{i},3));
end

p = num2str(task{2}{1}.randVars.p);
save(['stim_sltaskLoc' p],'mytex','scaledbI255')

%set new Gamma table for better luminance resolution sets the gamma table
%so that we can have finest possible control over the stimulus luminance
%thus contrast.
function setGammaTableForMaxContrast(maxContrast)
global stimulus;

%set reserved colors RGB triplets at the bottom of the table (first rows)
gammaTable(1:size(stimulus.colors.resvdTriplet,1),1:size(stimulus.colors.resvdTriplet,2)) = stimulus.colors.resvdTriplet;

%set gamma table
if maxContrast > 0
    
    %fill up the rest of the gamma table
    cmax = 0.5 + maxContrast/2;
    cmin = 0.5 - maxContrast/2;
    luminanceVals = cmin:((cmax - cmin)/(stimulus.colors.nPossCol-1)):cmax;
    
    %now get the linearized range The gamma table has 256 possible
    %luminance values from 0:1/255:1 We interpolate linearly the npossible
    %red/green/blue values at luminance values = luminanceVals
    redLinearized   = interp1(0:1/255:1,stimulus.linearizedGammaTable.redTable,luminanceVals,'linear');
    greenLinearized = interp1(0:1/255:1,stimulus.linearizedGammaTable.greenTable,luminanceVals,'linear');
    blueLinearized  = interp1(0:1/255:1,stimulus.linearizedGammaTable.blueTable,luminanceVals,'linear');
    
    %get rid of NaN
    if isnan(redLinearized(1))
        redLinearized(1) = redLinearized(2) - (redLinearized(3)-redLinearized(2));
        greenLinearized(1) = greenLinearized(2) - (greenLinearized(3)-greenLinearized(2));
        blueLinearized(1) = blueLinearized(2) - (blueLinearized(3)-blueLinearized(2));
    end
    
    %add these values to the table
    gammaTable((stimulus.colors.minStimIx:stimulus.colors.mxSIx+1),:) = [redLinearized;greenLinearized;blueLinearized]';
end
gammaTable(stimulus.colors.minStimIx,:)=[1 0 0];

%set gamma table
mglSetGammaTable(gammaTable);

%remember what current max contrast is that we can display
stimulus.currentMaxContrast = maxContrast;

%convert from mouse's to aperture's coord
function mouseinfo = pos2ap(fixStimulus,mouseinfo,stimulus)

%and the mouse coordinates relative to the center of the target screen
mouseinfo.x.pix = mouseinfo.x.pix-fixStimulus.pospix(1); %(in pix)
mouseinfo.y.pix = mouseinfo.y.pix-fixStimulus.pospix(2);

%convert pixels 2 angles (mglLines works with angles)
mouseinfo.x.angle.screen = mouseinfo.x.pix*mglGetParam('xPixelsToDevice');
mouseinfo.y.angle.screen = mouseinfo.y.pix*mglGetParam('yPixelsToDevice');

%calculate the transformation parameter that position the mouse on the
%aperture calculate length of the arrow (pythagoras theorem)
arrowsz2 = mouseinfo.x.angle.screen^2 + mouseinfo.y.angle.screen^2;%(angle)

%calculate transformation parameter
transpara = stimulus.radiusdeg^2/arrowsz2;

%transform actual coordinates to on-aperture coordinates.
mouseinfo.x.angle.onap=sqrt(mouseinfo.x.angle.screen^2*transpara)*sign(mouseinfo.x.angle.screen);
mouseinfo.y.angle.onap=sqrt(mouseinfo.y.angle.screen^2*transpara)*sign(mouseinfo.y.angle.screen);

%convert from PowerMate's to aperture's coord
function mouseinfo = PowerMate2ap(mouseinfo,stimulus,task)
global Inicoord;
%Set speed, initial position and size of arrow. wheelspeed=0.008;
wheelspeed=0.032;
%fixed initAngle=(225+270)/2*pi/180; random
initAngle.rad = task.thistrial.initAngledeg*pi/180;
r = stimulus.radiusdeg;

%Calculate arrow's coordinates on aperture
mouseinfo.x.angle.onap= r * cos(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap)); %(-) --> goes left when turn left, vvs.
mouseinfo.y.angle.onap= r * sin(initAngle.rad + wheelspeed*(mouseinfo.x.pix - Inicoord.x.pix.onap));

%convert from polar to cartesian coord
function coord = polar2cartesian(theta,r)
%theta is an angle in degree r is the radius of the unit circle Coord are
%in visual angle Record angle in degree
theta2.deg=theta;
%Convert from degree to radian
theta2.rad=theta2.deg*pi/180;
%Calculate visual angles coordinates
coord.x=r*cos(theta2.rad);
coord.y=r*sin(theta2.rad);
