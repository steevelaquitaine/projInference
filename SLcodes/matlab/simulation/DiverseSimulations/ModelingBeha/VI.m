
% Variational inference algorithm
%     author: Steeve Laquitaine
%       date: 13/06/05

%       goal: fit data with a Gaussian mixture model with two components
%         (bimodal distribution)
%             Our goal is to infer the posterior distribution for the means
%             covariance matrices and mixture coefficients, given a data set X.

%    keyword: clustering, classification, pattern recognition

%  reference: Pattern recognition and machine learning 2006, Bishop

%       note: VI maximizes the posterior (MAP estimate) distribution over
%       the parameters by weighting the likelihood computed by EM with a
%       prior over the parameters. This prevents overfitting.


% Bugs:
% ------
% Wes calculation is wrong.
% Calculation of Ses was wrong.
% Ses calculated for display is not the relevant one.
%


function VI( ~ )
%EM Summary of this function goes here
%   Detailed explanation goes here

% Set parameters of distribution
% number of data
n = 1000;
k = 2;
m = [[3 1.5] ; [0.3 0.8]];
s(:,:,1) = [0.8 0.3; 0.3 0.8];
s(:,:,2) = [1.2 0.0; 0.0 0.2];
p = [0.3 0.7];

% Draw from the mixture
[X, Z] = makeMixGauss(m, s, p, n, k);

% Plot
plotMixGauss(X, Z, 'mixgaussactual');
plotMeanStdEllipse(m,s,p);

% Variational inference
[mes ses pes] = gaussMixVI(X, 6);



% Make Variational inference
function [mes Ses pes] = gaussMixVI(X, k, dispFit)

% check default arguments
if nargin < 3,dispFit = 1;end
%keyboard

% get dimensionality of data X
D = size(X, 2);

% get sample size of data X
N =  length(X);

% Fixed hyperparameters: initial guess about the parameters of the
% priors over the model parameters.
% (see "pattern recognition and machine learning p102").
% (see http://en.wikipedia.org/wiki/Variational_Bayesian_methods)

% Dirichlet prior over mixture coefficients
a0 = ones(1,k)*(10^(-2.5));
pes = ones(1,k)*1/k;

% Gaussian-Wishart prior over means and covariances of the Gaussian
% components. The fixed initial hyperparameters are:
% B0: covariance matrix scaling parameter
% m0: means
% v0: is the number of degrees of freedom of the distribution
% w0: covariance matrices: is a D x D scale matrix.
B0 = ones(1,k);
m0 = rand(D,k);
v0 = 3;%D+1; %mathworks
w0(:,:,1) = eye(D);
w0(:,:,2) = eye(D);
w0(:,:,3) = eye(D);
w0(:,:,4) = eye(D);
w0(:,:,5) = eye(D);
w0(:,:,6) = eye(D);

% Initialize hyperparameters
aes = a0;
Bes = B0;
mes = m0;
ves = repmat(v0,1,k);
Wes = w0;
i =  1:D;

% bypass 0/0 case
small = 1e-10;

% VI's loop
StopCriterion = 0;
while StopCriterion ~= 1
    % StopCriterion = StopCriterion + 1
    
    % E-step: use the current variational distributions over the model
    % parameters to evaluate three moments and hence evaluate responsibilities.
    % xbar is 0 when responsibilities are 0.
    % -------------------------------------------------
    % expectations from Dirichlet
    Elogpi = psi(0, aes) - psi(0, sum(aes));%checked
    
    % expectations from Wishart
    stuff1 = D*log(2) + sum(psi(0,bsxfun(@minus,repmat(ves,2,1)+1,i')*0.5));
    for ki = 1 : k
        ElogLambda(ki) = stuff1(ki) + log(det(Wes(:,:,ki)));%checked
    end
    
    % get reponsibilities
    for ki = 1:k
        for ni = 1:N
            diff = X(ni,:)-mes(:,ki)';
            logRho(ni,ki) = Elogpi(ki) + 0.5*ElogLambda(ki) - D/(2*Bes(ki)) - (ves(ki)/2)*diff*Wes(:,:,ki)*diff';
        end
    end

    % normalize
    rnk = exp(logRho);
    rnk = rnk ./ repmat(sum(rnk,2),1,k);
    
    % samples of the sources
    Nk = sum(rnk) + small;
    for ki = 1 : k
        % new mean
        xbar(:,ki) = sum(bsxfun(@times, X, rnk(:,ki)),1) / Nk(ki);
        
        % and weighted covariances
        diffFromMean = (X-repmat(xbar(:,ki)',N,1))';
        Sk(:,:,ki) = (repmat(rnk(:,ki)',2,1).*diffFromMean)*diffFromMean'/Nk(ki);
    end
    
    
    % Variational M - step: update hyperparameters with responsibilities
    % -------------------------------------------------------------------
    Bes = B0 + Nk;
    mes = bsxfun(@rdivide, bsxfun( @plus, bsxfun(@times, Nk, xbar), bsxfun(@times, B0,m0)), Bes);
    for ki = 1 : k
        Wes(:,:,ki) = inv(inv(w0(:,:,ki)) + Nk(ki)*Sk(:,:,ki)+((B0(ki)*Nk(ki))/(B0(ki)+Nk(ki)))*(xbar(:,ki)-m0(:,ki))'*(xbar(:,ki)-m0(:,ki)));
    end
    ves = v0 + Nk;
    aes = a0 + Nk;
    
    
    % Get what we want for display
    % -----------------------
    for ki = 1:k
        Ses(:,:,ki) = inv(Wes(:,:,ki))/(ves(ki)-D-1); %Where did Justin found this equation ?
    end
    pes = Nk'/N;%see Hagai Attias' paper
    
    
    % Calculate variational lower bound
    % -----------------------
    % 1
    for ki  = 1 : k
        mahalkindof2(ki) = (xbar(:,ki) - mes(:,ki))'*Wes(:,:,ki)*(xbar(:,ki) - mes(:,ki));
        ElogpiXZuSk(ki) = Nk(ki)*(ElogLambda(ki) - D/Bes(ki) - ves(ki)*trace( Wes(:,:,ki)*Sk(:,:,ki)) - ves(ki)*mahalkindof2(ki) - D*log(2*pi));
    end
    ElogpiXZuS = 0.5*sum(ElogpiXZuSk);
    
    % 2
    ElogpiZp = sum(sum(bsxfun(@times, rnk, Elogpi)));
    
    % 3
    C = gamma(sum(a0))/prod(gamma(a0)); % checked
    Elogpip = log(C) + (a0(1)-1)*sum(Elogpi); % replace a0 by a unique a0;
    
    % 4
    PoG = prod((v0+1-i)/2);
    stuff2 = (v0 - D - 1)*sum(ElogLambda)/2;
    for ki  = 1 : k
        mahalkindof3(ki) = (mes(:,ki)-m0(:,ki))'*Wes(:,:,ki)*(mes(:,ki)-m0(:,ki));
        stuff3(ki) = D*log(B0(ki)/2*pi) + ElogLambda(ki) - D*B0(ki)/Bes(ki) - B0(ki)*ves(ki)*mahalkindof3(ki);
        stuff4(ki) = ves(ki)*trace(inv(w0(:,:,ki))*Wes(:,:,ki));
    end
    B = det(w0(:,:,1)).^(-v0/2) * (  2^(v0*D/2) * pi^(D*(D-1)/4) * PoG  )^-1; % replace a0 by a unique a0;
    stuff3 = sum(stuff3)/2;
    stuff4 = sum(stuff4)/2;
    ElogpiuS = stuff3 + k*log(B) + stuff2 - stuff4;
    
    % 5
    ElogqZ = sum(sum(rnk.*log(rnk))); %checked
    
    % 6
    Elogqp = sum((aes - 1).*Elogpi) + log(C);% checked
    
    % 7
    for ki=1:k
        PoG2(ki) = prod((ves(ki)+1-i)/2);
        B2(ki)= det(Wes(:,:,ki)).^(-ves(ki)/2) * (  2^(ves(ki)*D/2) * pi^(D*(D-1)/4) * PoG2(ki)  )^-1;
        stuff2b(ki) = (ves(ki) - D - 1)*ElogLambda(ki)/2;
    end
    HqS = -log(B2) - stuff2b + ves*D/2;
    ElogquS = sum((ElogLambda + D*log(Bes/2*pi) - D)/2 - HqS);
    
    % Lower bound
    LowerBound = ElogpiXZuS + ElogpiZp + Elogpip + ElogpiuS - ElogqZ - Elogqp - ElogquS;
    
    % plot
    % ----------------------------------------------
    % plot fit
    % --------------------------------
    % hard clustering
    % [dummy indicator] = max(rnk);
    % soft clustering
    indicator = rnk;
    
    plotMixGauss(X, indicator);
    titlestr = plotMeanStdEllipse(mes', Ses, pes);
    
    % title
    %     title(sprintf('%s\nlog like: %f', titlestr));%, logllh));
    drawnow
end
% draw nicely (with points color scaled by responsibility)
if dispFit
    disp(sprintf('(gaussmix) redisplaying with responsibilities'));
    % plot data
    cla;
    plotMixGauss(X, rnk(1,:));
    % plot fit
    titlestr = plotMeanStdEllipse(mes', Ses, pes);
    % title
    %title(sprintf('%s\nlog like: %f',titlestr,loglike));
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
plot(X(:,1), X(:,2), 'k.');

% draw point labeled according to which distribution they come from
subplot(1,2,2);
kn = max(Z');
% if kn is greater than one this is an indicator variable which
% says precisely which group each point is in, so color appropriately
if kn > 1.01
    for k = 1 : kn
        thisDist = find(Z==k);
        plot(X(thisDist,1), X(thisDist,2),'.','Color', getSmoothColor(k,kn,'hsv'));
        hold on
    end
else
    % otherwise a probability, so each points color gets scaled by the probability
    n = size(X,1);
    k = size(Z,2);
    % color1 = getSmoothColor(1,2,'hsv');
    % color2 = getSmoothColor(2,2,'hsv');
    %     for ni = 1:n
    %         %plot(X(ni,1), X(ni,2), '.', 'Color', color1*Z(ni) + color2*(1-Z(ni)));
    %         hold on
    %     end
    for ki = 1 : k
        colori(ki, :) = getSmoothColor(ki, k, 'hsv');
    end
    for ni = 1 : n
        colorthis = sum(bsxfun(@times, colori, Z(ni,:)'),1);
        plot(X(ni,1), X(ni,2), '.', 'Color', colorthis/max(colorthis));
        hold on
    end
end

% Plot ellipses
function titlestr = plotMeanStdEllipse(m, s, p)

k = length(s);
titlestr =[];

% Plot means
for ki = 1 : k
    plot(m(ki,1), m(ki,2), ...
        's', ...
        'MarkerFaceColor', getSmoothColor(ki,k,'hsv'), ...
        'MarkerEdgeColor','w',...
        'MarkerSize',16);
    
    % Plot std ellipses
    ellipse(s(2,2,ki), -s(1,2,ki), s(1,1,ki), 1, m(ki,:), getSmoothColor(ki,k,'hsv'),12)
end

m=m';
if numel(p) == 6
    titlestr = sprintf('m1: [%0.2f %0.2f] m2: [%0.2f %0.2f] m3: [%0.2f %0.2f] m4: [%0.2f %0.2f] m5: [%0.2f %0.2f] m6: [%0.2f %0.2f] \n s1: [%0.2f %0.2f] s2: [%0.2f %0.2f] s3: [%0.2f %0.2f] s4: [%0.2f %0.2f] s5: [%0.2f %0.2f] s6: [%0.2f %0.2f] \n p1: %0.2f p2: %0.2f p3: %0.2f p4: %0.2f p5: %0.2f p6: %0.2f \n',...
        m(1,1),m(2,1), m(1,2),m(2,2), m(1,3),m(2,3), m(1,4),m(2,4), m(1,5),m(2,5), m(1,6),m(2,6), ...
        s(1,1,1),s(2,2,1), s(1,1,2),s(2,2,2), s(1,1,3),s(2,2,3), s(1,1,4),s(2,2,4), s(1,1,5),s(2,2,5), s(1,1,6),s(2,2,6),...
        p(1), p(2), p(3), p(4), p(5), p(6));
end

title(titlestr);

% bypass matlabi numerical overfolow
function s = logsumexp(x, dim)
% Compute log(sum(exp(x),dim)) while avoiding numerical underflow.
%   By default dim = 1 (columns).
% Written by Michael Chen (sth4nth@gmail.com).
if nargin == 1,
    % Determine which dimension sum will use
    dim = find(size(x)~=1,1);
    if isempty(dim), dim = 1; end
end

% subtract the largest in each column
y = max(x,[],dim);
x = bsxfun(@minus,x,y);
s = y + log(sum(exp(x),dim));
i = find(~isfinite(y));
if ~isempty(i)
    s(i) = y(i);
end