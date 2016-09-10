
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
%load data03
%load data04
%load dataAll

load data01_control


%Input k
%coherence:{'0.06','0.12','0.06'};
%prior:{'80','40','20','10'};
k(1:3)=[27.7 5 0.7];
k(4:7)=[0.7 1.5 13.6 219];
k(8)=0.05;
k(9)=14.5;
k(10)=1;
k(11)=mean([0.000043 0.27 1.45 1.55]);

