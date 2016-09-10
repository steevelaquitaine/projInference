
% SLmakeSound.m
%
%     author: steeve laquitaine
%       date: 140625
%    purpose: make sound

function SLmakeSound(nSeconds)

%Samples per second
Fs = 1000;      
 
%Tone frequency, in Hertz
toneFreq = 300; 
y=sin(linspace(0, nSeconds*toneFreq*2*pi, round(nSeconds*Fs)));

%sound
sound(y, Fs); 