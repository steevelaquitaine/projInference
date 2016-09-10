

%SLmakeTuningCurves.m
%
%       $Id: SLmakeTuningCurves.m $
%        by: steeve laquitaine
%      date: 141106
%   purpose: Create N neurons tuning neurons with different feature
%            preference (e.g., motion direction)
%
%     usage:
%
%           o.nNeuron = 100;
%           o.maxResp = 20;
%           o.baseFR = 1;
%           o.kTuning = 4;
%           o.TuningDiSpace = 1 : 1 : 360;
%
%           o = SLmakeTuningCurves(o)
%
%Description
%
%   - N neurons with tuning responses over motion directions
%
%   - Tuning curves are von Mises (f can be a non integer because 
%     it is an average firing rate over trials)
%     fi(direction) =  Gain*von mises + base firing;
%
%   - Neurons have preferred directions between 0 and 360º.
%
%
%reference: 
%     http://homepages.inf.ed.ac.uk/pseries/CCN14/lab4.pdf

function o = SLmakeTuningCurves(o)

%feature preferences
o.pref = round(360/o.nNeuron : 360/o.nNeuron : 360);

%feature space (e.g., direction)
o.D = length(o.TuningDiSpace);

%Tuning
o.f = nan(o.D,o.nNeuron);

%vectorize for speed
TunDi = SLde2r(o.TuningDiSpace,1)'; 
TunDiAll = TunDi(:,ones(1,o.nNeuron)); %direction space (rad)

pref = SLde2r(o.pref,1);
prefAll = pref(ones(o.D,1),:); %neurons preferred directions (rad)

o.f = o.maxResp * exp(o.kTuning*cos(TunDiAll - prefAll) - o.kTuning) + o.baseFR; %tuning


