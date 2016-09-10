
%Initialize model parameters for Bayesian inference with learnt and
%cardinal prior
tic
close all; 
clear
load data01

%Input k
%k(1:3):llh strengths
%k(4:7):learnt prior strength
%k(8)  :cardinal prior strengths
%k(9)  :probability to choose randomly
%k(10) :motor noise
%k(11) :probability to switch to learnt prior (=0 for no switch)
clear k
k(1:3)=[80 40 20];
k(4:7)=[1.74 4.77 10.74 34.25];
k(8)=0;
k(9)=0.001;
k(10)=15;
k(11)=0;
