
%SLsimPPCafterAdaptation.m
%
%       $Id: SLsimPPCafterAdaptation.m $
%        by: steeve laquitaine
%      date: 141108
%
%   purpose: Create a Probabilistic Population Code. 
%            Motion directions (e.g., a space of feature) are encoded in 
%            neural tuning responses and likelihood over motion directions 
%            to the readout of an actually presented motion direction.
%            Tuning responses are adapted.
%           
%     usage: write the function before calling a figure
%
%          SLshadesOfGreyLinspecer


function C = SLshadesOfGreyLinspecer

a = log(1:0.01:10)';

%scale between 0 and 1
a = (a-min(a)) / max(a);

C = [a a a];
axes('NextPlot','replacechildren', 'ColorOrder',C);
