%author: Steeve Laquitaine
  %date: 140331 last modification 140422
 %usage: [estimatedDir,Ctest,CtrEachDir]=reconstructDir(Btest,W,Ctr,motionDirTr)   
         %BtestL: matrix of Bold response (m voxels, n trials)
              %W: matrix of channels weights (m voxels, k channels)
            %Ctr: matrix of channels outputs (k channels, n trials)
    %motionDirTr: vector of displayed motion directions

%Description: decode or reconstruct displayed motion directions from voxels 
%pattern of BOLD responses, channels weights and a lookup table of channel 
%outputs for each displayed direction used to train the forward encoding 
%model.

%Note 140407:
%- More voxels reduces variability of matrix of weights' conditionability.
%tested with 1000 neurons/voxels, does not seems to reduce conditionability
%itself.

%- More neurons in a voxel ..... matrix of weights' conditionability.
%tested with 5000 voxels.


function [estimatedDir,Ctest,CtrEachDir,Condition]=reconstructDir(Btest,W,Ctr,motionDirTr)
%matrix of channels outputs for test data (such that Btest=W*Ctest)
%W'*W must be invertible to have a unique solution. Otherwise infinity of
%solution because of multicollinearity.i.e., columns of X are linearly
%dependent.
%W'*W must be full rank i.e, rank(WmkLIP'*WmkLIP)=size(WmkLIP'*WmkLIP)
%det(WmkLIP'*WmkLIP) must be far from 0
%if matrix is full rank it must be well-conditioned: cond(W'*W) must be 
%close to 1.
%Matlab decides like that: cond(Wtmp'*Wtmp)<1/((max(size(Wtmp'*Wtmp)))*eps)

%Ctest=inv(W'*W)*W'*Btest;
if rank(W'*W)~=size(W'*W)
    fprintf('\n %12s \n','(reconstructDir) IMPORTANT WARNING: The matrix of weights is not full rank, thus not invertible. The output channels solution are not unique')
else
    fprintf('\n %12s \n','(reconstructDir) The matrix of weights is invertible but may be ill-conditioned. Make sure output variable "Condition" ~ 1')
end
Ctest=(W'*W)\W'*Btest;
Condition=cond(W'*W);

%known channel outputs for each motion direction from trainingsize()
[eachDir,IA]=unique(motionDirTr,'first');
CtrEachDir=Ctr(:,IA);
numDirTrained=size(CtrEachDir,2);

%compare estimated channel outputs "Ctest"  in each trial to known channel 
%outputs for each trained direction "CtrEachDir"
%Comparison is done for each trial for each trial and estimated direction
%the direction that produces max correlation between channel outputs.
numDirTested=size(Ctest,2);
parfor i=1:numDirTested
    Rbkp=[];
    for j=1:1:numDirTrained
        [R,P]=corrcoef(Ctest(:,i),CtrEachDir(:,j));
        Rbkp=[Rbkp R(2,1)];
    end
    [maxRbkp,pos]=max(Rbkp);
    maxR(i)=maxRbkp;
    estimatedDir(i)=eachDir(pos);
end
fprintf('\n %12s \n','(reconstructDir) reconstructing parameter...done')

