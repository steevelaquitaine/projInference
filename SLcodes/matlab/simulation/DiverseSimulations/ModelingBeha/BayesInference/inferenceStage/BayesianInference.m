


function BayesInference
%BAYESIANINFERENCE Summary of this function goes here
%   Detailed explanation goes here

% Create distributions


%---------------------------------------------------------
%% Strong motion - high certainty in the prior
%---------------------------------------------------------

subplot(2,2,1)
axis square
hold all

% highly certain (or informed) prior (gaussian)
% parameters of the prior
prior.mean=140;
prior.std=sqrt(1000);
obs.nb=10000000;
% sampledir=0:360;
sampledir=-360:360;


prior.values = prior.mean + prior.std.*randn(obs.nb,1);
prior.countperdir=hist(prior.values,sampledir); % check if mistake here
%------------------------------------------------------------------------
% % REPRESENTATIONS
% ------------------------------------------------------------------------
% evaluate probabilities
prior.probaperdir=prior.countperdir/obs.nb; 

% Check if probabilities sum to 1
sum(prior.probaperdir)

% Likelihood (gaussian)
% Set parameters of the Likelihood
% mean
evidence.mean=90;
% std
evidence.std=sqrt(300);
evidence.values = evidence.mean + evidence.std.*randn(obs.nb,1);
evidence.countperdir=hist(evidence.values,sampledir);

% Evaluate probabilities
evidence.probaperdir=evidence.countperdir/obs.nb;

% Check if probabilities sum to 1
sum(evidence.probaperdir)

% ------------------------------------------------------------------------
% % OPERATIONS
% ------------------------------------------------------------------------
% posterior
% p(b/a)=p(a/b)*p(a)/p(b)
% p(b)=p(b/a=1) + p(b/a=2) +....
posterior.probaperdir=prior.probaperdir.*evidence.probaperdir/sum(prior.probaperdir.*evidence.probaperdir);

% Check if probabilities sum to 1
sum(posterior.probaperdir)

% MARKERS
% Plot prior's best estimate (mean)
plot([prior.mean prior.mean],[0 0.04],'--b',...
    'linewidth',1)
% Plot evidence's best estimate (mean)
plot([evidence.mean evidence.mean],[0 0.04],'--k',...
    'linewidth',1)

% Plot distributions
plot(sampledir,posterior.probaperdir,'color',[.5 .5 .5],...
    'linewidth',2)
plot(sampledir,prior.probaperdir,'color','r',...
    'linewidth',2)
plot(sampledir,evidence.probaperdir,'color','k',...
    'linewidth',2)
ylim([0 0.04])
xlim([-100 300]) %replace 3 by 0 here after correcting mistake

%---------------------------------------------------------
%% Strong motion - high uncertainty in the prior
%---------------------------------------------------------
subplot(2,2,3)
axis square
hold all

% prior (gaussian)
% parameters of the prior
prior.mean=140;
prior.std=sqrt(4000);
obs.nb=1000000;
% sampledir=-360:360;

prior.values = prior.mean + prior.std.*randn(obs.nb,1);
prior.countperdir=hist(prior.values,sampledir);

% evaluate probabilities
prior.probaperdir=prior.countperdir/obs.nb;

% check if equal 1
sum(prior.probaperdir)

% likelihood (gaussian)
% parameters of the evidence
evidence.mean=90;
evidence.std=sqrt(300);

evidence.values = evidence.mean + evidence.std.*randn(obs.nb,1);
evidence.countperdir=hist(evidence.values,sampledir);

% evaluate probabilities
evidence.probaperdir=evidence.countperdir/obs.nb;

% check if equal 1
sum(evidence.probaperdir)


% posterior
% p(b/a)=p(a/b)*p(a)/p(b)
% p(b)=p(b/a=1) + p(b/a=2) +....
posterior.probaperdir=prior.probaperdir.*evidence.probaperdir/sum(prior.probaperdir.*evidence.probaperdir);

% check if equal 1
sum(posterior.probaperdir)

% MARKERS
% plot prior's best estimate (mean)
plot([prior.mean prior.mean],[0 0.04],'--b',...
    'linewidth',1)
% plot evidence's best estimate (mean)
plot([evidence.mean evidence.mean],[0 0.04],'--k',...
    'linewidth',1)


% plot distributions
plot(sampledir,posterior.probaperdir,'color',[.5 .5 .5],...
    'linewidth',2)
plot(sampledir,prior.probaperdir,'color','r',...
    'linewidth',2)
plot(sampledir,evidence.probaperdir,'color','k',...
    'linewidth',2)
ylim([0 0.04])
xlim([-100 300]) %replace 3 by 0 here after correcting mistake


%---------------------------------------------------------
%% Weak motion - high certainty in the evidence
%---------------------------------------------------------
subplot(2,2,2)
axis square
hold all

% prior (gaussian)
% parameters of the prior
prior.mean=140;
prior.std=sqrt(1000);
obs.nb=1000000;
% sampledir=-360:360;

prior.values = prior.mean + prior.std.*randn(obs.nb,1);
prior.countperdir=hist(prior.values,sampledir);

% evaluate probabilities
prior.probaperdir=prior.countperdir/obs.nb;

% check if equal 1
sum(prior.probaperdir)

% likelihood (gaussian)
% parameters of the evidence
evidence.mean=90;
evidence.std=sqrt(1000);

evidence.values = evidence.mean + evidence.std.*randn(obs.nb,1);
evidence.countperdir=hist(evidence.values,sampledir);

% evaluate probabilities
evidence.probaperdir=evidence.countperdir/obs.nb;

% check if equal 1
sum(evidence.probaperdir)

% posterior
% p(b/a)=p(a/b)*p(a)/p(b)
% p(b)=p(b/a=1) + p(b/a=2) +....
posterior.probaperdir=prior.probaperdir.*evidence.probaperdir/sum(prior.probaperdir.*evidence.probaperdir);

% check if equal 1
sum(posterior.probaperdir)

% MARKERS
% plot prior's best estimate (mean)
plot([prior.mean prior.mean],[0 0.04],'--b',...
    'linewidth',1)
% plot evidence's best estimate (mean)
plot([evidence.mean evidence.mean],[0 0.04],'--k',...
    'linewidth',1)


% plot distributions
plot(sampledir,posterior.probaperdir,'color',[.5 .5 .5],...
    'linewidth',2)
plot(sampledir,prior.probaperdir,'color','r',...
    'linewidth',2)
plot(sampledir,evidence.probaperdir,'color','k',...
    'linewidth',2)
ylim([0 0.04])
xlim([-100 300]) %replace 3 by 0 here after correcting mistake
%---------------------------------------------------------
%% Weak motion - high uncertainty in the evidence
%---------------------------------------------------------
subplot(2,2,4)
axis square
hold all

% prior (gaussian)
% parameters of the prior
prior.mean=140;
prior.std=sqrt(4000);
obs.nb=1000000;
% sampledir=-360:360;
% 
prior.values = prior.mean + prior.std.*randn(obs.nb,1);
prior.countperdir=hist(prior.values,sampledir);

% evaluate probabilities
prior.probaperdir=prior.countperdir/obs.nb;

% check if equal 1
sum(prior.probaperdir)

% likelihood (gaussian)
% parameters of the prior
evidence.mean=90;
evidence.std=sqrt(1000);
obs.nb=1000000;

evidence.values = evidence.mean + evidence.std.*randn(obs.nb,1);
evidence.countperdir=hist(evidence.values,sampledir);

% evaluate probabilities
evidence.probaperdir=evidence.countperdir/obs.nb;

% check if equal 1
sum(evidence.probaperdir)

% posterior
% p(b/a)=p(a/b)*p(a)/p(b)
% p(b)=p(b/a=1) + p(b/a=2) +....
posterior.probaperdir=prior.probaperdir.*evidence.probaperdir/sum(prior.probaperdir.*evidence.probaperdir);

% check if equal 1
sum(posterior.probaperdir)

% MARKERS
% plot prior's best estimate (mean)
plot([prior.mean prior.mean],[0 0.04],'--b',...
    'linewidth',1)
% plot evidence's best estimate (mean)
plot([evidence.mean evidence.mean],[0 0.04],'--k',...
    'linewidth',1)


% plot distributions
plot(sampledir,posterior.probaperdir,'color',[.5 .5 .5],...
    'linewidth',2)
plot(sampledir,prior.probaperdir,'color','r',...
    'linewidth',2)
plot(sampledir,evidence.probaperdir,'color','k',...
    'linewidth',2)
ylim([0 0.04])
xlim([-100 300]) %replace 3 by 0 here after correcting mistake


end

