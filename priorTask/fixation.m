%fixation.m
%
%       $Id: fixStairInitTask.m 708 2010-03-01 03:38:19Z justin $
%     usage: [fixTask myscreen]=fixStairInitTask(myscreen)
%        by: justin gardner
%      date: 09/07/06
% copyright: (c) 2006 Justin Gardner (GPL see mgl/COPYING)
%   purpose: Implements a fixation task. In this task, the fixation cross
%            starts out cyan and then darkens twice. The fixation cross then
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
%Modified by A Meso May 2012 for a fixation synched to stimuli but without
%a related task to solve

%Modified by Steeve Laquitaine August 2012 for a fixation synched to
%stimuli.

%Notes:
%Fixation point is a disk of 0.2 degrees radius in visual angle.

function [task myscreen]=fixation(myscreen)

%check arguments
if ~any(nargin == [1])
  help fixation
  return
end

%create the stimulus for the experiment, use defaults if they are
%not already set
global fixStimulus;
myscreen=initStimulus('fixStimulus',myscreen);

%fixation parameters
fixStimulus.fixWidth=3; %(cm) at this size, fixation also acts as a cardinal axis reference... 
fixStimulus.fixLineWidth=100; %degree
fixStimulus.pos=[0 0]; %positon in the center of the screen

%init the task
%Before the 131220 this seglen was set at inf; but it nows produces an
%error. It may be due to mgl functions update. I set it now at a large
%number it seems to work fine.
task{1}.seglen=100000;
[task{1} myscreen]=initTask(task{1},myscreen,@fixStartSegmentCallback,@fixDrawStimulusCallback);

%function that gets called at the start of each segment
function [task myscreen]=fixStartSegmentCallback(task, myscreen)

global fixStimulus;
%if this is the inter stimulus interval
  %training mode, clear screen here
%  if fixStimulus.trainingMode,mglClearScreen;end 
  fixStimulus.thisColor=127; %grey so that it is not too strong over many trials

  
%function that gets called every frame udpate to draw the fixation 
function [task myscreen]=fixDrawStimulusCallback(task, myscreen)
global fixStimulus;

%Stimulus size (degree of visual angle)
mglGluDisk(0,0,0.2,0.5,24,2)