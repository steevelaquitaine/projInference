


%slinitRunDual2AFC.m


%vertical orientG1 * orientG2 * respCue * 

function slinitRunDual2AFC

%grating angles
%(Baruni, from 3 deg to 45 from category boundary (45 deg))
p.orientG1 = [42:-3:0 48:3:90]; %bottom left
p.orientG2 = [42:-3:0 48:3:90]; %top right

p.orientG1(p.orientG1==45) = [];
p.orientG2(p.orientG1==45) = [];

%grating contrast
p.contrast = 1;

%response cue location
%1:top 0: bottom
p.respCue = [1 2]: