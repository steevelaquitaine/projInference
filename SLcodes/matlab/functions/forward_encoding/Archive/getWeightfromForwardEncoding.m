
 %author: steeve laquitaine
   %date: 140326
%purpose: get the forward encoding channels' weights for each voxel
%e.g., if I have 6 channels and 100 voxels' response I should have a matrix
%of 6 rows * 100 columns.


%input should be:
%Ctr (kxn): matrix of predicted channel outputs (value) for each channel 
%(row) and each trial (n)

%Btr (mxn) is matrix of BOLD response amplitude (value) for each voxel (row)
%for each trial (col).

%m: voxel number
%n: trials number
%k: channels number


function W=getWeightfromForwardEncoding(Btr,Ctr)

%Linear regression with matrix algebra produces the weights
%W relates the matrix of BOLD responses with the Channels
%Basically Btr(theta)=W.*Ctr(theta)
W=Btr*(Ctr')*((Ctr*Ctr')^-1);
