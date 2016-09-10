% gaussmix.m
%
%        $Id:$ 
%      usage: gaussmix()
%         by: justin gardner
%       date: 05/29/13
%    purpose: 
%
function retval = testmixofgauss()

% check arguments
if ~any(nargin == [0])
  help gaussmix
  return
end

% parameters of distribution
n = 1000;
k = 2;
m = [[0 1.5] ; [0.3 0.8]];
s(:,:,1) = [0.8 0.3;0.3 0.4];
s(:,:,2) = [1.2 0.0;0.0 0.2];
p = [0.3 0.7];

% get draws from this this mixture distribution
[X Z] = makeMixGauss(n,k,m,s,p);

% plot them
plotMixGauss(X,Z,'mixgaussactual');
plotMeanStdEllipse(m,s);

% estimate the means, std and mixing coefficients using EM
[mest sest pest] = gaussMixEM(X,2);


%%%%%%%%%%%%%%%%%%%%
%    gaussMixEM    %
%%%%%%%%%%%%%%%%%%%%
function [means s p] = gaussMixEM(X,k,dispFit)

% check default arguments
if nargin < 3,dispFit = 1;end

keyboard

stopCriteria = 0;
while stopCriteria ~= 1

  % display fit
  if dispFit
    % plot data
    [dummy indicator] = max(rnk');
    cla;plotMixGauss(X,indicator); 
    % plot fit
    titlestr = plotMeanStdEllipse(means,s);
    % title
    title(sprintf('%s\nlog like: %f',titlestr,loglike));
    drawnow
  end
end

% draw nicely (with points color scaled by responsibility)
if dispFit
  disp(sprintf('(gaussmix) redisplaying with responsibilities'));
  % plot data
  cla;plotMixGauss(X, rnk(:,1));
  % plot fit
  titlestr = plotMeanStdEllipse(means,s);
  % title
  title(sprintf('%s\nlog like: %f',titlestr,loglike));
  drawnow
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    plotMeanStdEllipse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titlestr = plotMeanStdEllipse(means,s)

k = length(s);
% plot means
for ki = 1:k
  plot(means(ki,1),means(ki,2),'s','MarkerFaceColor',getSmoothColor(ki,k,'hsv'),'MarkerEdgeColor','w','MarkerSize',16);
  % plot std ellipses
  ellipse(s(2,2,ki),-s(1,2,ki),s(1,1,ki),1,means(ki,:),getSmoothColor(ki,k,'hsv'),12)
end

titlestr = sprintf('m1: [%0.2f %0.2f] m2: [%0.2f %0.2f] s1: [%0.2f %0.2f] s1 cov: %0.2f s2: [%0.2f %0.2f] s2 cov: %0.2f',means(1,1),means(1,2),means(2,1),means(2,2),s(1,1,1),s(2,2,1),s(1,2,1),s(1,1,2),s(2,2,2),s(2,1,2));
title(titlestr);

%%%%%%%%%%%%%%%%%%%%%%
%    plotMixGauss    %
%%%%%%%%%%%%%%%%%%%%%%
function plotMixGauss(X,Z,figname)

if nargin < 3
  mlrSmartfig('gaussmix','reuse');
  clf;
else
  mlrSmartfig(figname,'reuse');
  clf;
end

% draw points all mixed together
subplot(1,2,1);
plot(X(:,1),X(:,2),'k.');

% draw point labeled according to which distribution they come from
subplot(1,2,2);
kn = max(Z);
% if kn is greater than one this is an indicator variable which
% says precisely which group each point is in, so color appropriately
if kn>1
  for k = 1:kn
    thisDist = find(Z==k);
    plot(X(thisDist,1),X(thisDist,2),'.','Color',getSmoothColor(k,kn,'hsv'));
    hold on
  end
else
  % otherwise a probability, so each points color gets scaled by the probability
  color1 = getSmoothColor(1,2,'hsv');
  color2 = getSmoothColor(2,2,'hsv');
  n = size(X,1);
  for ni = 1:n
    plot(X(ni,1),X(ni,2),'.','Color',color1*Z(ni)+color2*(1-Z(ni)));
    hold on
  end
end
%%%%%%%%%%%%%%%%%%%%%%
%    makeMixGauss    %
%%%%%%%%%%%%%%%%%%%%%%
function [X Z] = makeMixGauss(n,k,m,s,p)

% this function returns n random sampled values from a mixture of k 2D gaussian distributions
% with means specified as a (k,2) arrays of means in the two dimensions
% and s is a cell of length k with elements that are the (2 x 2) covariance matrices
% p (k x 1) is the mixing coefficients of the distributions (i.e. what proportion of points come
% from what distribution)

% which distribution each point will come from, choose to match the mixing coefficients
% normalize p to make sure it is a probability
p = p / sum(p);
% take cumulative sum since this will make it easy to randomly select values
p = cumsum([0 p]);

% now draw randomly so that Z is a number from 1-k for each of the n samples (this is using 
% the terminology for the latent variable Z which is an indicator that specifies which dist
% each point is from)
r = rand(1,n);
for ki = 1:k
  Z(r >= p(ki)) = ki;
end

X = nan(n,2);
for ki = 1:k
  % get how many points in this distribution
  ni = sum(Z==ki);
  
  % now compute values from a 2D gaussian with the appropriate
  % mean and covariance
  R = chol(squeeze(s(:,:,ki)));
  X(Z==ki,:) = repmat(m(ki,:),ni,1) + randn(ni,2)*R;
end