% gaussmix.m
%
%        $Id:$ 
%      usage: gaussmix()
%         by: justin gardner
%       date: 05/29/13
%    purpose: Test program for fiting a mixture of gaussian to data
%             First generates a mixture of gaussians data set
%             Then fits with an EM algorithm or variational Bayes 
%             Based on Bishop book
%
function retval = gaussmix()

% check arguments
if ~any(nargin == [0])
  help gaussmix
  return
end

% parameters of distribution
n = 1000;
k = 2;
m = [[0 1.5] ; [0.3 0.8]];
%m = [[0 1.5] ; [2.8 0.2]];
s(:,:,1) = [0.8 0.3;0.3 0.4];
s(:,:,2) = [1.2 0.0;0.0 0.2];
p = [0.3 0.7];

% get draws from this this mixture distribution
[X Z] = makeMixGauss(n,k,m,s,p);

% plot them
plotMixGauss(X,Z,'mixgaussactual');
plotMeanStdEllipse(m,s);

% estimate the means, std and mixing coefficients using EM
%[mest sest pest] = gaussMixEM(X,2);

% estimate the means, std and mixing coefficients using variational bayes
% this also returns a structure with full parameters of the posterior function
% estimated by variational bayes
[mest sest pest params] = gaussMixVariationalBayes(X,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    gaussMixVariationalBayes    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [means s p params] = gaussMixVariationalBayes(X,k,dispFit)

% check default arguments
if nargin < 3,dispFit = 1;end

% number of dimensions in X
n = size(X,1);
D = size(X,2);

%%%%%%%%%%%%%%%%%%%%
% Hyper parameters. 
%%%%%%%%%%%%%%%%%%%
% These are all parameters of the prior distribution for the quantities that we are trying
% to estimate - i.e. the prior distribution for the Z latent variable (which distribution
% each point is taken from) , and the mean and cov of the gaussian distributions from which 
% each point is modeled to be chosen from.

% alpha is for the Dirichlet distribution used as a prior for the categorical variable Z
% which specifies which distribution each data point comes from. It is a variable of length
% k to specify an alpha value for each distribution. It roughly corresponds to the "number of
% observations associated with each component of the mixture. Making this small will make the posteriro
% primarily influenced by the data rather than prio (p475). p480 discusses what happens as you change 
% this value larger. Basically you favor more components.
alpha0 = ones(1,k)*(10^(-2.5));

% m0 and beta0 are used for specifing the normal distribution from which is the prior for the mean of
% each distribution. We set the mean of these priors to all be 0 and beta0 to ?
%[indicator m0] = getKmeans(X,k,dispFit);
m0 = rand(k,D);
%m0 = zeros(k,D);
beta0 = ones(1,k);

% W0 and v0 are used for the Wishart distribution which specifies the prior for the covariance
% of the gaussian distributions. What should these be set to?
W0 = eye(D);
v0 = 3;
invW0 = inv(W0);

% set the current parameters to the initial parameters
alphak = alpha0;
betak = beta0;
mk = m0;
for ki = 1:k
  Wk(:,:,ki) = W0;
end
vk = repmat(v0,k,1);

stopCriteria = 0;
while stopCriteria ~= 1
  
  % E-step: Calculate the responsibilities - i.e. the probability each data point comes from
  % each of the distribution. Analgous to E-step in EM. 

  % Equation 10.66 log pi tilde which specifies the expeceted log mixing
  % coefficients according to the digamma distribution
  logPiTilde = psi(0,alphak)-psi(0,sum(alphak));
%  disp(sprintf('Mixing: %s',mynum2str(exp(logPiTilde))));
  
  % Equation 10.65 expected value of the log of the determinant of the precision matrix
  for ki = 1:k
    logLambdaTilde(ki) = sum(psi(0,(vk(ki)+1-(1:D))/2))+D*log(2)+log(det(Wk(:,:,ki)));
  end

  % Equation 10.67. Note that this is computing the responsibilities (the probability
  % that each point comes from each of the gaussians. Note the samiliarity to EM where
  % this essentially is computing the probabilty (according to the gaussina distribution)
  % that each point comes from that distribution. There are other factors in here because
  % you don't just have the best (max likelihood) estimate in each step, but the full posterior
  % distribution, so you have to integrate over all possible parameter values to compute
  % this probabilty. This is equation 10.67. Note that we calculate the log of each
  % side of the equation 10.67 this is to avoid some numerical issues with trying
  % to exponentiate to a very big negative values
  for ki = 1:k
    for ni = 1:n
      diff = X(ni,:)-mk(ki,:);
      logRho(ni,ki) = logPiTilde(ki) + 0.5*logLambdaTilde(ki) - D/(2*betak(ki)) - (vk(ki)/2)*diff*Wk(:,:,ki)*diff';
    end
  end
  % normalize the responsabilities to make them into probabilities
  rnk = exp(logRho);
  rnk = rnk ./ repmat(sum(rnk,2),1,k);

  % M-step: Now that we have the responsiblites calculate the best parameter estimates - but
  % this does not mean to compute just the best (maximum likelihood) parameters, but the 
  % full likelihood function which is parameterized by a gauss-wishart distribution. That is 
  % the expected value of the gauss is the best mean and the expected value of the wishart
  % is the best standard-deviation of the gaussian. So, we are in this case optimizing the
  % various parameters of the gaussian-wishart in this step

  % this is the effective N for this k (i.e. the sum of the responsiblities for this distribution);
  for ki = 1:k
    % compute Nk and add small fundge term to avoid 0 Nk values
    Nk(ki) = sum(rnk(:,ki)) + 1e-10;

    % and calculate weighted mean
    xMeank(ki,:) = sum(repmat(rnk(:,ki),1,D).*X)/Nk(ki);
    
    % and weighted covariance
    diffFromMean = (X-repmat(xMeank(ki,:),n,1))';
    s(:,:,ki) = (repmat(rnk(:,ki)',2,1).*diffFromMean)*diffFromMean'/Nk(ki);
  end

  % update parameters of distributions specifying parameters
  alphak = alpha0 + Nk;
  betak = beta0+Nk;
  vk = v0 + Nk;
  for ki = 1:k
    % mk is the mean of the gauss distribution of means. Note that this is
    % just a weighted average of the prior mean (m0) and the current mean (xMeank)
    % since betak is just beta0+Nk as defined above.
    mk(ki,:) = (1/betak(ki)).*(beta0(ki).*m0(ki,:) + Nk(ki)*xMeank(ki,:));
    % calculate the new Wk
    Wk(:,:,ki) = inv(invW0 + Nk(ki)*s(:,:,ki)+((beta0(ki)*Nk(ki))/(beta0(ki)+Nk(ki)))*(xMeank(ki,:)-m0(ki,:))'*(xMeank(ki,:)-m0(ki,:)));
    % calculate the best estimte for the standard deviation (for displaying)
    sest(:,:,ki) = inv(Wk(:,:,ki))/(vk(ki)-D-1);
  end

  % display fit
  if dispFit
    % plot data
    plotMixGauss(X,rnk,'gaussmixvarbayes');
    % plot fit
    titlestr = plotMeanStdEllipse(mk,sest,exp(logPiTilde));
    xaxis(-4,6);yaxis(-1,3);
    % title
    title(sprintf('Variational Bayes\n%s',titlestr));
    drawnow
  end
end


%%%%%%%%%%%%%%%%%%%%
%    gaussMixEM    %
%%%%%%%%%%%%%%%%%%%%
function [xMeank s pi] = gaussMixEM(X,k,dispFit)

% check default arguments
if nargin < 3,dispFit = 1;end

% number of dimensions in X
n = size(X,1);
D = size(X,2);

% first initialize estimates with kmeans (note that we could use the matlab built in kmeans
% function for this, but its more fun to implement ourselves ...
[indicator xMeank] = getKmeans(X,k,dispFit);
%[indicator xMeank] = kmeans(X,k);

% initialize standard deviations with the standard deviation across all the data
sall = cov(X);
for ki = 1:k
  s(:,:,ki) = sall;
end

% ok, now cycle through
lastloglike = -inf;
epsilon = 10e-8;
stopCriteria = 0;
while stopCriteria ~= 1
  % E step - calculate what the expectation of the responsibilites are, i.e. which distribution
  % each point is likely to have come for
  for ki=1:k
    rnk(:,ki) = mvnpdf(X,xMeank(ki,:),s(:,:,ki));
  end
  % scale to make into a probability
  rnk = rnk./repmat(sum(rnk,2),1,k);
  
  % M step - calculate the means given the responsiblities calculated above
  % first calculate means
  for ki = 1:k
    % compute Nk and add small fundge term to avoid 0 Nk values
    Nk(ki) = sum(rnk(:,ki)) + 1e-10;

    % and calculate weighted mean
    xMeank(ki,:) = sum(repmat(rnk(:,ki),1,D).*X)/Nk(ki);
    
    % and weighted covariance
    diffFromMean = (X-repmat(xMeank(ki,:),n,1))';
    s(:,:,ki) = (repmat(rnk(:,ki)',2,1).*diffFromMean)*diffFromMean'/Nk(ki);
  end

  % calculate log likelihood of fit
  for ki = 1:k
    p(ki,:) = (Nk(ki)/n).*mvnpdf(X,xMeank(ki,:),s(:,:,ki));
  end
  % get log likelihood, note that we sum over the k dimension
  % since the probablities are or (i.e. the same point could be
  % in the 1st or the 2nd or the kth distribution) and we multiply
  % over the data points since that is an and. This is implemented
  % as two sums (the second happens after taking a log so is a prod)
  loglike = sum(log(sum(p)));

  % check stopping criteria
  stopCriteria = (loglike-lastloglike) < epsilon;
  lastloglike = loglike;

  % display fit
  if dispFit
    % plot data
    plotMixGauss(X,rnk,'gaussmixem');
    % plot fit
    titlestr = plotMeanStdEllipse(xMeank,s);
    % title
    title(sprintf('EM: %s\nlog like: %f',titlestr,loglike));
    drawnow
  end
end

% draw nicely (with points color scaled by responsibility)
if dispFit
  disp(sprintf('(gaussmix) redisplaying with responsibilities'));
  % plot data
  cla;plotMixGauss(X,rnk,'gaussmixem',1);
  % plot fit
  titlestr = plotMeanStdEllipse(xMeank,s);
  % title
  title(sprintf('EM: %s\nlog like: %f',titlestr,loglike));
  drawnow
end

% return estimated mixing coefficients 
pi = Nk/n;


%%%%%%%%%%%%%%%%%%%
%    getKmeans    %
%%%%%%%%%%%%%%%%%%%
function [indicator means] = getKmeans(X,k,dispFit)

% check default arguments
if nargin < 3,dispFit = 0;end

% size of data vector
n = size(X,1);
indicator = ones(n,1);

% initialize indicator to random
indicator = ceil(rand(n,1)*k);

stopCriteria = 0;
lastMeans = nan;
while stopCriteria ~= 1

  % Do M-step first - calculate means
  for ki = 1:k
    means(ki,:) = mean(X(indicator==ki,:));
  end

  % Do E-Step: recalculate indicator as whichever mean is nearest to each point
  dist = [];
  for ki = 1:k
    dist(ki,:) = sum((X-repmat(means(ki,:),n,1)).^2,2);
  end
  [dummy indicator] = min(dist);

  % display
  if dispFit
    % display the distributions
    plotMixGauss(X,ind2oneofk(indicator,k),'gaussmixkmeans');
    % display the means
    for ki = 1:k
      plot(means(ki,1),means(ki,2),'o','MarkerFaceColor',getSmoothColor(ki,k,'hsv'),'MarkerEdgeColor','w','MarkerSize',16);
    end
    title(sprintf('kMeans: m1=[%0.2f %0.2f] m2=[%0.2f %0.2f]',means(1,1),means(1,2),means(2,1),means(2,2)));
    drawnow
  end
  
  % set stopping criteria, as when the means no longer change
  if isequal(lastMeans,means),stopCriteria = 1;end
  lastMeans = means;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    plotMeanStdEllipse    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titlestr = plotMeanStdEllipse(means,s,pi)

% deal with input arguments.
% titleDist tells which two distributions to display info about in title
titleDist = [1 2];
if nargin>=3
  [dummy maxN] = sort(pi);
  titleDist = [maxN(end) maxN(end-1)];
  dispMixing = 1;
else
  pi = ones(1,size(means,1));
  dispMixing = 0;
end
k = length(s);
epsilon = 10e-8;
% plot means
for ki = 1:k
  % only plot if mixing coefficent greater than a small number
  if pi(ki) > epsilon
    % plot means
    plot(means(ki,1),means(ki,2),'s','MarkerFaceColor',getSmoothColor(ki,k,'hsv'),'MarkerEdgeColor','w','MarkerSize',16);
    % plot std ellipses
    ellipse(s(2,2,ki),-s(1,2,ki),s(1,1,ki),1,means(ki,:),getSmoothColor(ki,k,'hsv'),12)
  end
end

% get the means and distributions for which
% we wish to display information (usually the ones
% with the largest effective n
mean1 = means(titleDist(1),:);
mean2 = means(titleDist(2),:);
std1 = s(:,:,titleDist(1));
std2 = s(:,:,titleDist(2));

% display information
titlestr = sprintf('m1: [%0.2f %0.2f] m2: [%0.2f %0.2f] s1: [%0.2f %0.2f] s1 cov: %0.2f s2: [%0.2f %0.2f] s2 cov: %0.2f',mean1(1),mean1(2),mean2(1),mean2(2),std1(1,1),std1(2,2),std(1,2),std2(1,1),std2(2,2),std2(2,1));
% display the mixing coefficients
if dispMixing
  titlestr = sprintf('%s\npi = %s',titlestr,mynum2str(pi));
end
% dispaly title
title(titlestr);

%%%%%%%%%%%%%%%%%%%%%%
%    plotMixGauss    %
%%%%%%%%%%%%%%%%%%%%%%
function plotMixGauss(X,Z,figname,blendColors)

if (nargin < 3) || isempty(figname)
  mlrSmartfig('gaussmix','reuse');
  clf;
else
  mlrSmartfig(figname,'reuse');
  clf;
end
if nargin < 4, blendColors = 0;end

% draw points all mixed together
subplot(1,2,1);
plot(X(:,1),X(:,2),'k.');

% draw point labeled according to which distribution they come from
subplot(1,2,2);
kn = size(Z,2);
% if kn is greater than one this is an indicator variable which
% says precisely which group each point is in, so color appropriately
if ~blendColors
  [dummy maxZ] = max(Z,[],2);
  for k = 1:kn
    thisDist = find(maxZ==k);
    plot(X(thisDist,1),X(thisDist,2),'.','Color',getSmoothColor(k,kn,'hsv'));
    hold on
  end
else
  % otherwise a probability, so each points color gets scaled by the probability
  for ki = 1:kn
    color{ki} = getSmoothColor(ki,kn,'hsv');
  end
  n = size(X,1);
  for ni = 1:n
    thisColor = [0 0 0];
    for ki = 1:kn
      thisColor = thisColor + Z(ni,ki)*color{ki};
    end
    plot(X(ni,1),X(ni,2),'.','Color',thisColor);
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
  Zindexed(r >= p(ki)) = ki;
end

X = nan(n,2);
for ki = 1:k
  % get how many points in this distribution
  ni = sum(Zindexed==ki);
  
  % now compute values from a 2D gaussian with the appropriate
  % mean and covariance
  R = chol(squeeze(s(:,:,ki)));
  X(Zindexed == ki,:) = repmat(m(ki,:),ni,1) + randn(ni,2)*R;
end
Z = ind2oneofk(Zindexed,k);

%%%%%%%%%%%%%%%%%%%%
%    ind2oneofk    %
%%%%%%%%%%%%%%%%%%%%
function Zoneofk = ind2oneofk(Zindexed,k)

% get number of points
n = length(Zindexed(:));

% convert from an index to a 1-of-k representation where there are k 
% columns all filled with 0 except one column - that column number represents
% which distribution that point comes from
Z = zeros(n,k)';
Z(Zindexed+(0:k:(k*(n-1)+k-1))) = 1;
Zoneofk = Z';

%%%%%%%%%%%%%%%%%%%%
%    oneofk2ind    %
%%%%%%%%%%%%%%%%%%%%
function Zindexed = oneofk2ind(Zoneofk)

% convert from one-of-k to an indexed

% get number of rows
n = size(Zoneofk,1);
% find one values
[row ind] = find(Zoneofk);
% init the zindexed
Zindexed = zeros(1,n);
% and set every row to the column that had a one
Zindexed(row) = ind;
