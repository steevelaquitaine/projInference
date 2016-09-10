
%Simulate forward encoding (FE)
%Simulate BOLD responses from a Bayesian neural network (training 
%and test data), run forward encoding model on training data and use the 
%channels weights to reconstruct motion directions in test data.

%Observations:
%- Simulation shows that decoding gets better with more voxels (1000's good)
%- Increasing trial repetition reduces variability of decoded direction.
%- priorstd must be 0.0727 to have same strength as likelihood.

%1- Get weights of FE channels from simulated BOLD data MTBOLDmn (m voxels, 
%n trials) to train the FE model. Be sure to have repetitions otherwise the
%weights will be poorly trained.
close  all
%motionDirTrain=repmat(1:1:360,1,7);
motionDirTrain=repmat(1:1:360,1,1);
[MTBOLDmnTrain,LIPBOLDmnTrain]=simulateBOLD(motionDirTrain,...
    'priorstd=0',...
    'numVoxels=100',...
    'numVoxelsLIP=100',...
    'NeurInaVoxelMT=100000',...
    'NeurInaVoxelLIP=100000',...
    15,...
    'displayNetwork=off',...
    'displayBOLD=on');

%2- run forward encoding model
%get a matrix of channels' output Ctrkn and the weights of each channels
%for each voxel.
[CtrknMT,channelsMT,dirSpaceMT,preferredDirMT,WmkMT]=doForwardEncoding(MTBOLDmnTrain,6,motionDirTrain,'display=off');
[CtrknLIP,channelsLIP,dirSpaceLIP,preferredDirLIP,WmkLIP]=doForwardEncoding(LIPBOLDmnTrain,6,motionDirTrain,'display=off');

%3- decode motion direction used to train FE 
%from MT data
estimatedDirMT=reconstructDir(MTBOLDmnTrain,WmkMT,CtrknMT,motionDirTrain);
figure; hold all; drawCircStat(estimatedDirMT',motionDirTrain');

%from LIP data
estimatedDirLIP=reconstructDir(LIPBOLDmnTrain,WmkLIP,CtrknLIP,motionDirTrain);
figure; hold all; drawCircStat(estimatedDirLIP',motionDirTrain');

% %5- simulate test data
% % motionDirTest=repmat(1:1:360,1,7);
% motionDirTest=repmat(1:1:360,1,1);
% [MTBOLDmnTest,LIPBOLDmnTest]=simulateBOLD(motionDirTest,...
%     'priorstd=0.0727',... 'numVoxels=1000',... 'numVoxelsLIP=1000',15,...
%     'display=off');
% 
% %6- decode motion direction from LIP data
% estimatedDirLIPtest=reconstructDir(LIPBOLDmnTest,WmkLIP,CtrknLIP,motionDirTrain);
% figure; hold all;
% drawCircStat(estimatedDirLIPtest',motionDirTest');