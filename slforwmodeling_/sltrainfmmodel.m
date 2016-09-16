
%sltrainfmmodel.m
%
%author : steeve laquitaine


function fm = sltrainfmmodel(b,svec)

%get channel responses from 5 direction tuned channels
%average channel responses to vector of stimulus directions
%Ni instances x 8 channels
fm = slsimfmchannels;
f_k_s = fm.f_k_s';
C = fm.f_k_s(:,svec)';

%Ordinary linear regression
Wtrained = pinv(C)*b;
fm.Wtrained = Wtrained';
