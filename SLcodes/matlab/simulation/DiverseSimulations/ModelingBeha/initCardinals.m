
%Initialize model parameters
tic
close all; 
clear
load data01

%Input k
%k(1:3):llh strengths
%k(4)  :cardinal prior strength
%k(5)  :probability to choose randomly
%k(6)  :motor noise
%k(7)  :probablity to switch to prior mean (must be 0)
clear k
k(1:3)=[167 16 4];
k(4)=0;
k(5)=0.0;
k(6)=67.2;
k(7)=0;
