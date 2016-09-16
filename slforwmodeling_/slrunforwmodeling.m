

%slrunforwmodeling.m
%
%
%author : steeve laquitaine
%  date : 160915
%purpose: model brain voxel responses with a forward model that assumes
%         voxel responses are the linear sum of 5 sinewave "channel" 
%         function with 5 different direction preference (phase) raised
%         to power 4.
%
%
%
%simulate training instance matrix "b_train" : 
%626 instances x 320 voxels matrix of responses to 5 different stimulus 
%directions where 
clear;
Ni = 626;
Nv = 320;
s = [15 85 155 225 295];
svec = randsample(s,Ni,'true');
W = rand(Nv,5);
[b_train,C] = slsimfminstances(svec,W);

%% train the model's weights
fm = sltrainfmmodel(b_train,svec);
slfmPlotWeights(fm.Wtrained,W,fm.phi_k);

%% decode from test dataset and get average channel responses
%re-center the channel responses relative to the match between the channel 
%preference and the displayed orientation and average them
svec_test = randsample(s,Ni,'true');
b_test = slsimfminstances(svec_test,W);
fm = sltestfmmodel(b_test,svec_test,fm);
%plot
figure('color','w')
plot(s-fm.center,fm.meantestCresp,'.-','markersize',20)
box off



