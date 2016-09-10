function [ val_flct ] = DecayGssRndWlk(m,values,iteration, lambda,theta,sigma )
% Steeve laquitaine 
% started and ended 13/06/2011

%DECAYGSSRNDWLK Summary of this function goes here
% Set a decaying gaussian random walk
% Bibliography: Daw et al, Nature Neuroscience, 2006; Wunderlich et al,
% PNAS, 2011 

% DETAILED EXPLANATION GOES HERE
% m: number of values (stimuli)
% values (what are the initial values)
% iteration (number of trials)
% lambda (decay parameter )
% theta (decay center)
% sigma (for v the zero-mean Gaussian diffusion noise)

lambda=0.9836; % Wunderlich et al;  Daw et al (0.8)
theta=500; % random rewarding
sig=2.8*10; % corrupting noise's parameter

val_flct(1,:)=values;
for t=1:iteration-1
    v(t,:)=randn(m,1)*sig'; % zero-mean gaussian
    val_flct(t+1,:)=lambda * val_flct(t,:) + (1 - lambda) * theta +  v(t,:);
end

end

