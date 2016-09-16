

%slrunforwmodelingKFoldCV.m
%
%
%author : steeve laquitaine
%  date : 160915
%purpose: model brain voxel responses with a K fold cross-validated 
%         forward model that assumes
%         voxel responses are the linear sum of 5 sinewave "channel" 
%         function with 5 different direction preference (phase) raised
%         to power 4.
%
%
%
%simulate training instance matrix "b_train" : 
%626 instances x 320 voxels matrix of responses to 5 different stimulus 
%directions where 
Ni = 626;
Nv = 320;
s = [15 85 155 225 295];
svec = randsample(s,Ni,'true');
W = rand(Nv,5);
[instances,C] = slsimfminstances(svec,W);

%% cross-validated decoding of a dataset
fm = slvoxfmKFoldCVdec(instances,svec,5);