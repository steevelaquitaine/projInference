
%SLmakePopulationResponses.m
%
% author: steeve laquitaine
%   date: 150631
%purpose: Simulate N neurons responses (spike counts) to a stimulus feature (e.g., motion direction).
%
% usage :
%
%       o.nNeuron = 40;         %number of different tunings
%       o.rSpace = 0:1:60       %spike counts space
%       o.PRgivenDi             %a length(rSpace) *. s Stim *. N neurons matrix of probabilities
%       o.stim = 45;
%       o = SLmakePopulationResponses(o,Optiondisplay)

function o = SLmakePopulationResponses(o,Optiondisplay)

%check
if o.stim == 0
    fprintf('(SLmakePopulationResponses) o.stim = 0. o.stim should be defined on 1:1:360')
    keyboard
end

%get neurons response to stimulus: sample a spike count at stimulus
%position on tuning function
o.rs = nan(1,o.nNeuron);

for Ni = 1 : o.nNeuron
    
    o.rs(Ni) = randsample(o.rSpace,1,'true',o.PRgivenDi(:,o.stim,Ni));
    
end

%display
if nargin >1
    o = display2(Optiondisplay,o);
end