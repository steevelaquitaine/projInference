
%sltrainPPmodel.m


function Wtrained = sltrainPPmodel(b,C)

%Ordinary linear regression
Wtrained = pinv(C)*b;
Wtrained = Wtrained';
