
%Initialize model parameters

%The fraction of trials k(8) for random estimation and k(10) for response 
%bias input to simulate data are different from the one implemented by 
%the simulation because trials are assigned to processes stochastically.
%The real fractions of trials are "kreal".

%to simulate full Bayes inference set k(10)=0;

%initialize settings for subjects 3
tic
close all; 
clear
load data01
%load data02
%load data03
%load data04
%load dataAll

%Input k
k(1:3)=[39.8 7.4 2.9];
k(4:7)=[0.09 2.3 2.2 4];
k(8)=0.0007;
k(9)=67.2;
k(10)=0;