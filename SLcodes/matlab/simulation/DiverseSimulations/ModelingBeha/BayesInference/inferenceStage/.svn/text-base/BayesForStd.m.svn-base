function [ output_args ] = BayesForStd( input_args )

%% What this does:
% Compute the qualitative predictions made about the posterior std by
% bayesian inference


% Steps
% - Combine two gaussians (prior and likelihood)
% - Normalize the resulting posterior to probabilities.
% - Calculate the std of the posterior
% - Increase the distance between the prior's mean and the likelihood mean
% - How does the std behave as the distance increases?


% states (directions displayed)
x = 5: 1: 355;


% prior
P.u = 225; %mean
P.s = 20;  %std
P.p = 1./(P.s*sqrt(2*pi)) * exp(-0.5*((x - P.u)/P.s).^2)'; %density
P.p = P.p/sum(P.p);%probability


% likelihood
l.u = 5; 
l.s = 20;
l.p = 1./(l.s*sqrt(2*pi)) * exp(-0.5*((x - l.u)/l.s).^2)';
l.p = l.p/sum(l.p);


% Draw
figure('color','w')
hold on
plot(P.p,'w')
area(x,P.p,'FaceColor',[.6 .6 .6],'edgecolor','k');

plot(l.p,'w')
area(x,l.p,'FaceColor',[.2 .2 .2],'edgecolor','k');
xlim([min(x) max(x)])











