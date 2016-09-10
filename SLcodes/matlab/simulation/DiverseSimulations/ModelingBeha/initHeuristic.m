
%Initialize model parameters

%The fraction of trials k(8) for random estimation and k(10) for response 
%bias input to simulate data are different from the one implemented by 
%the simulation because trials are assigned to processes stochastically.
%The real fractions of trials are "kreal".

%Description:
%-ratio kp/kl controls bias.
%-amplitude of kp &kl controls variability.
%-fraction of trials with random estimation, k(8)
%-motor noise, k(9)
%-Beta (k(10)) increases the probability that the stronger representation 
%win the competition (prior or llh). It increases the steepness of softmax
%sigmoid function. 
% Beta=0 is unbiased

%initialize settings for subjects 3
tic
%close all; 
clear
%load data01
%load data02
load data03HeuristinitControl.mat
%load data04
%load dataAll

%Input k
%coherence:{'0.06','0.12','0.06'};
%prior:{'80','40','20','10'};
k(1:3)=[5299 14.7 2.6];
k(4:7)=[0.8 2.7 14.5 160];
k(8)=0.007;
k(9)=30;
k(10)=1;
