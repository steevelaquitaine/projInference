%fixationLoc
%
%author: steeve laquitaine

function [task,myscreen] = fixationLoc(myscreen)


%check arguments
if ~any(nargin == [1])
  help fixation
  return
end

%create the stimulus for the experiment, use defaults if they are
%not already set
global fixStimulus;
fixStimulus.pospix = [myscreen.screenWidth myscreen.screenHeight]/2;
myscreen = initStimulus('fixStimulus',myscreen);

%fixation parameters
fixStimulus.fixWidth=3; %(cm) at this size, fixation also acts as a cardinal axis reference... 
fixStimulus.fixLineWidth=100; %degree
fixStimulus.pos=[0 0]; %positon in the center of the screen

%init the task
%Before the 131220 this seglen was set at inf; but it nows produces an
%error. It may be due to mgl functions update. I set it now at a large
%number it seems to work fine.
task{1}.seglen = 100000;
[task{1} myscreen] = initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback);

%function that gets called at the start of each segment
function [task, myscreen] = fixStartSegmentCallback(task, myscreen)
global stimulus
global fixStimulus;
%if this is the inter stimulus interval
  %training mode, clear screen here
%  if fixStimulus.trainingMode,mglClearScreen;end 
fixStimulus.thisColor = myscreen.fixCr;

  
%function that gets called every frame udpate to draw the fixation 
function [task myscreen]=fixDrawStimulusCallback(task, myscreen)
global stimulus
global fixStimulus;

%Stimulus size (degree of visual angle)
mglGluDisk(0,0,0.2,myscreen.fixCr,24,2)