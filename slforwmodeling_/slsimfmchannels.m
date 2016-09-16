

%slsimfmchannels.m
%
%
%author : steeve laquitaine
%purpose: simulate 5 channels (half rectified cosines 
%         raised to power 4) over direction circular space 
%         1:1:360 degs with direction preferences
%         [15 85 155 225 295].
%
%usage :
%
%       pp = slsimfmchannels
%
%
%
%outputs :
%
%           pp.f_k_s : matrix of 8 channels tuning functions
%                      (K channels x S stimulus directions)
%           pp.phi_k : channels direction preferences
%               pp.K : number of channels    


function fm = slsimfmchannels

%stimulus direction space
fm.s = 1:1:360;

%assumed 8 tuning functions
%with at least one at 180 deg
fm.phi_k = [15 85 155 225 295]';
fm.K = length(fm.phi_k);

%directions tuning functions 
%s - phi_k for each phi_k and s
phi_k = fm.phi_k;
s = fm.s;
sTophi_k = bsxfun(@(s,phi_k) s-phi_k,s,phi_k);
fm.f_k_s = max( 0,cos(pi*(sTophi_k)/180) ).^5;
