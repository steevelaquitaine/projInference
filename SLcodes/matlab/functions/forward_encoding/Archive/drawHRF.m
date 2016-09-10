
%author: steeve Laquitaine
  %date: 140422
 %usage: drawHRF
%Description: draw hemodynamic response function

 %reference:https://labs.psych.ucsb.edu/ashby/gregory/Matlab.html

function drawHRF

%Clear all workspaces
clear all; close all; clc;
figure('color','w')

%Set the parameter values of the hrf
T0=0; 
n=4; 
lamda=2;
 
%time axis
t=0:.01:30;
 
%hrf
hrf=((t-T0).^(n-1)).*exp(-(t-T0)/lamda)/((lamda^n)*factorial(n-1));
 
%Plot
plot(t,hrf); axis([0 30 0 .12]);
box off