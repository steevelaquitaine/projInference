
%Initialize model parameters
tic
close all; 
clear
load data01

%Input k
%k(1:3):llh strengths
%k(4:7):learnt prior strengths
%k(8)  :probability to choose randomly
%k(9)  :motor noise
%k(10) :probability of flat likelihood
%k(11) :cardinal prior strength
clear k
k(1:3)=[0 0 0];
k(4:7)=[0.1 2 8 33];
k(8)=0.0007;
k(9)=67.2;
k(10)=0;
k(11)=0;
