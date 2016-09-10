
% EM algorithm
%     author: Steeve Laquitaine, inspired from Justin Gardner' code
%       date: 13/06/05

%       goal: fit data with a Gaussian mixture model with two components 
%         (bimodal distribution)
%    keyword: clustering, classification, pattern recognition
%  reference: Pattern recognition and machine learning 2006, Bishop

%       note: As expected the algorithm sometimes stays stuck. One
%       component takes responsibility for all the data and the other
%       component collapses. It is overfitting. 
%       EM maximizes log likelihood that the Gaussian mixture explains the
%       data.




function EM( input_args )
%EM Summary of this function goes here
%   Detailed explanation goes here


% n: number of samples
% k: number of components (gaussians)
% m: mean (2D)
% s: covariance matrix
% p: mixture coefficient


% Set parameters of distribution
n = 1000;
k = 2;
m = [[3 1.5] ; ...% compo 1
    [0.3 0.8]];     % compo 2

s(:,:,1) = ...      % compo 1
    [0.8 0.3; ... 
    0.3 0.4];     

s(:,:,2) = ...      % compo 2
    [1.2 0.0; ...
    0.0 0.2];

p = [0.3 0.7];


% Get draws from this mixture of Gaussians
[X, Z] = makeMixGauss(m, s, p, n, k);

% Draw
plotMixGauss(X, Z, 'mixgaussactual');
plotMeanStdEllipse(m,s);

% Expectation maximization
[mes ses pes] = gaussMixEM(X, 2);

% Fit means, std and mixing coefficients using EM
function [mes ses pes] = gaussMixEM(X, k)

% check default arguments
if nargin < 3,dispFit = 1;end

keyboard

% Step 1: Set some initial parameters
mes = [[0 3] ; ...
    [2 2]];       

ses(:,:,1) = ...      
    [0.5 0.3; ... 
    0.3 0.2];     

ses(:,:,2) = ...      
    [1.8 0.0; ...
    0.0 0.5];

pes = [0.5 0.5];

logllh = 10000;
while logllh ~= 0
    
    % E step: Calculate responsibilities with current parameters
    % Get likelihood of data
    % loop over component
    for ki = 1 : k;
        % get distributions
        Pxnk(ki,:) = make2DGaussPdf(X, mes(ki,:), ses(:, :, ki));
    end
    
    % get llh
    LLH = bsxfun(@times, Pxnk, pes');
    
    % get responsibilities.
    rnk = bsxfun(@rdivide, LLH, sum(LLH));
    
    % M step: Re-estimate the parameters using the current responsibilities
    N =  length(X);
    Nk = sum(rnk,2);
    
    for ki = 1 : k
        % new mean
        mes(ki,:) = sum(bsxfun(@times, X', rnk(ki,:)),2)/Nk(ki);
        
        % new covariance matrix
        distance = bsxfun(@minus, X, mes(ki,:));
        for n = 1 : N
            sesn(:,:,n,ki) = rnk(ki,n)*distance(n,:).'*distance(n,:);
        end
        ses(:,:,ki) = sum(sesn(:,:,:,ki), 3)/Nk(ki);
        
        % new mixture coefficients
        pes(:, ki) = Nk(ki)/N;
    end
    
    % Step 4: evaluate log likelihood
    % get distributions
    for ki = 1 : k;
        Pxnk(ki,:) = make2DGaussPdf(X, mes(ki,:), ses(:, :, ki));
    end
    
    % get LLH
    LLH = bsxfun(@times, Pxnk, pes');
    % get log LLH
    logllh = sum(log(sum(LLH))); 
    
    
    % plot fit
    % --------------------------------
    % hard clustering
%     [dummy indicator] = max(rnk);
    % soft clustering
    indicator = rnk(1,:);
    
    plotMixGauss(X, indicator);
    titlestr = plotMeanStdEllipse(mes, ses);
    
    % title
    title(sprintf('%s\nlog like: %f', titlestr, logllh));
    drawnow  
end
% draw nicely (with points color scaled by responsibility)
if dispFit
  disp(sprintf('(gaussmix) redisplaying with responsibilities'));
  % plot data
  cla;
  plotMixGauss(X, rnk(:,1));
  % plot fit
  titlestr = plotMeanStdEllipse(means,s);
  % title
  title(sprintf('%s\nlog like: %f',titlestr,loglike));
  drawnow
end

% Make Gaussian mixtures
function [X, Z] = makeMixGauss(m, s, p, n, k)

% n: number of samples
% k: number of components (gaussians)
% m: mean (2D)
% s: covariance matrix
% p: mixture coefficient


% Make sure p is a probability
p = p / sum(p);

% take cumulative sum to make it easy to randomly select values
% This maybe works only with bimodal gaussians
p = cumsum([0 p]);

% Make latent variable Z that switches on or off each sources k 
% (components). Typically Z is binary {0,1}. Value of z depends on mixture
% components p such that P(Z=1) = p.
% Here for simplicity Z directly indicates the components activated.
r = rand(1, n);
% loop over components
for ki = 1 : k
  Z(r >= p(ki)) = ki;
end

% Draw data from 2D Gaussian
X = nan(n,2);
% loop over components
for ki = 1 : k
    % get number of points produced by this component
    ni = sum(Z==ki);
    
    % now compute values from a 2D gaussian
    X(Z==ki,:) = mvnrnd(m(ki,:), s(:,:,ki), ni);
end    

% Plot Gaussian mixtures 
function plotMixGauss(X, Z, figname)

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
if kn > 1
    for k = 1 : kn
        thisDist = find(Z==k);
        plot(X(thisDist,1), X(thisDist,2),'.','Color',getSmoothColor(k,kn,'hsv'));
        hold on
    end
else
    % otherwise a probability, so each points color gets scaled by the probability
      color1 = getSmoothColor(1,2,'hsv');
      color2 = getSmoothColor(2,2,'hsv');
    n = size(X,1);
    for ni = 1:n
        plot(X(ni,1), X(ni,2), '.', 'Color', color1*Z(ni) + color2*(1-Z(ni)));
        hold on
    end
end

% Plot ellipses
function titlestr = plotMeanStdEllipse(means, s)

k = length(s);

% Plot means
for ki = 1 : k
    plot(means(ki,1), means(ki,2), ...
      's', ...
      'MarkerFaceColor', getSmoothColor(ki,k,'hsv'), ...
      'MarkerEdgeColor','w',...
      'MarkerSize',16);
  
  % Plot std ellipses
  ellipse(s(2,2,ki), -s(1,2,ki), s(1,1,ki), 1, means(ki,:), ...
      getSmoothColor(ki,k,'hsv'),12)
end

titlestr = sprintf('m1: [%0.2f %0.2f] m2: [%0.2f %0.2f] s1: [%0.2f %0.2f] s1 cov: %0.2f s2: [%0.2f %0.2f] s2 cov: %0.2f',means(1,1),means(1,2),means(2,1),means(2,2),s(1,1,1),s(2,2,1),s(1,2,1),s(1,1,2),s(2,2,2),s(2,1,2));
title(titlestr);

% Make 2D Gaussian
function Px = make2DGaussPdf(X, m, s)

% X: data (2D)
% m: mean (2D)
% s: covariance matrix(DxD matrix)

% Set two dimensions
D = 2;

% Get pdf
for i = 1 : length(X)
    % get mahalanobis distance (use inner product)
    mahalanobisDis(i) = (X(i,:)-m) * s^(-1) * (X(i,:)-m)';
    % get pdf
    Px(i) = ((2*pi)^(-D/2)) * (det(s)^(-0.5)) * exp(-0.5*mahalanobisDis(i));
end




