
%SLsearchWJMBayesinitP.m
%
% author: steeve laquitaine
%   date: 150710 last modif 150728
%purpose: get reasonable kappa and kp initial parameters for WJM model fit
%         plot logL for each parameter combination of the grid
%         10*10*10 kappas (for each coherence)
%         grid search over kappa initial parameters (no km)
%         search for two basics models: 
%           - 1) accurate priors (fixed true kappa_s)
%           - 2) flat priors (fixed kappa_s = 0)
%
% usage: 
%
%       initPset = SLmkWJMinitPSet(kappa24set,kappa12set,kappa06set);
%        o = SLsearchWJMBayesinitP({'sub01'},initPset,'wjmLogl201to1000fitpGrid.mat','experiment','vonMisesPrior');
%
%
%
%You can also select a set of the initial parameters
%
%       initPsetAll = SLmkWJMinitPSet(0:10:90,0:10:90,0:10:90);
%       initPset = SLchooseNinitPSet(initPsetAll,1:1000);
%       o = SLsearchWJMBayesinitP({'sub01'},initPset,'wjmLogl201to1000fitpGrid.mat','experiment','vonMisesPrior');


function o = SLsearchWJMBayesinitP(subjects,initPset,filename,varargin)

tic 

%backup .m file 
mfilename = SLgetActivemFile;
myMfile = SLbackupMfileInWSpace(mfilename);
o.mfilename  = mfilename;
o.myMfile = myMfile;

%use parallel processing for speed up
delete(gcp)   %shutdown current parallel pool
matlabpool 10 %run a new one

%-----------------
%set path and data
%-----------------
fprintf('%s \n','(SLfitBayesianModel_WJM) Gathering data ...')
fprintf('%s \n','(SLfitBayesianModel_WJM) Please set dataPath ...')
dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');
%dataPath = input('Enter project path (e.g., /dataPsychophy/proj01...) : ','s');%request path
cd(dataPath)
varg = [varargin,'dataPath',dataPath]; %experiment and path
db   = SLMakedatabank(subjects,varg);
data         = db.estimatesDeg;      %subjects estimates
pstd         = db.Pstd;              %task priors std
sStrg        = db.stimStrength;      %stimulus coherence
s            = db.stimFeatureDeg;    %stimulus direction
[taskC,~,posC] = SLuniqpair([pstd sStrg s]); %task conditions

%get initial parameter set
fitKappa = initPset;
nGrid = size(fitKappa,1);
    
%loop over grid combinations of initial parameters
fitpAll = [];
parfor i = 1 : nGrid
    
    fitp = [fitKappa(i,:) 0.748 2.77 8.7 33.3 1e-6 3e30]; %true prior
    %fitp = [fitKappa(i,:) 0 0 0 0 1e-6 3e30]; %flat prior
    
    pm0 = taskC;
    pm0(taskC==0.24) = fitp(1);
    pm0(taskC==0.12) = fitp(2);
    pm0(taskC==0.06) = fitp(3);
    pm0(taskC==80)   = fitp(4);
    pm0(taskC==40)   = fitp(5);
    pm0(taskC==20)   = fitp(6);
    pm0(taskC==10)   = fitp(7);
    prand            = fitp(8);
    km               = fitp(9);
    s       = SLde2r(taskC(:,3),0); %s
    c       = taskC(:,2);           %coh
    kappa   = pm0(:,2);             %sensory strength
    kappa_s = pm0(:,1);
    
    %get logL
    logl(i) = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fitp,posC);
    fitpAll = [fitpAll; fitp]
end

tocs = toc;
o.duration = tocs;
o.logl = logl;
o.fitp = fitpAll;
save(filename)

%get logl of data
function [loglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC)

%Simulating measurements x
Ndots = 328;      % number of dots
Ntrials = 5000;   % number of trials

%prior parameters
Ns = 100;
svec = linspace(1,360,Ns+1)';
svec = SLde2r(svec(1:end-1),0);  %vector of hypothesized motion directions

%predicted estimates space
predEst = 1:1:360;
nCon = max(posC);

%each condition
parfor Ci = 1 : nCon
        
    %measurements
    coherent = (rand(Ntrials,Ndots) < c(Ci)); %which dots are moving coherently
    x = circ_vmrnd(s(Ci)*ones(Ntrials,Ndots),kappa(Ci));
    x = coherent .* x + (1-coherent) .* rand(Ntrials,Ndots)*2*pi; % Ntrials by Ndots
    
    %prior
    Mean_s = SLde2r(225,0);          %known by subj (fixed)
    prior = exp(kappa_s(Ci) * cos(svec - Mean_s));
    
    %Likelihoods L(s) = p(x|s)
    x = permute(repmat(x, [1 1 Ns]),[3 1 2]); % to make dimensions match: Ns by Ntrials by Ndots
    like_per_dot = (1-c(Ci))/2/pi + c(Ci)/2/pi/besseli(0,kappa(Ci)) * exp(kappa(Ci) * cos(bsxfun(@minus,x,svec))); %(1-c)/2pi + c*V(svec,x,kappa)
    like = prod(like_per_dot,3);
    
    %Posteriors p(s|x)
    posterior = bsxfun(@times, prior, like); % Ns by Ntrials
    
    %Estimation: posterior means (circular)
    xtemp = sum(bsxfun(@times, cos(svec), posterior));
    ytemp = sum(bsxfun(@times, sin(svec), posterior));
    posteriormeans = atan2(ytemp,xtemp); % 1 by Ntrials)
    posteriormeans = SLra2d(posteriormeans);
    
    %estimate distributions
    PestGivModel(Ci,:) = hist(posteriormeans,predEst);
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    PestGivModel(Ci,:) = (1-prand)*PestGivModel(Ci,:) + prand*(1/360);               %add random estimation
    PestGivModel(Ci,:) = SLcircConv(PestGivModel(Ci,:),vmPdfs(0:1:359,0,km,'norm')');%add motor noise

    fprintf('%s %i %.2f %.2f %.2f  %.3f %s \n','Cond: [',Ci,s(Ci),c(Ci),kappa(Ci),kappa_s(Ci),']')
end

%%%need check that condition data match conditions predictions

%get logl
for Ci = 1 : nCon
    PestGivModel(Ci,PestGivModel(Ci,:)<=0) = 10^-320;                   %get rid of negatives due to convolution
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    
    %logl
    thisCon = posC==Ci;
    logl(Ci) = sum(log(PestGivModel(Ci,data(thisCon))));             
end
loglAll = sum(logl);
fitP = fp; %best fitp


%sample circular density
function Y = circ_vmrnd(theta, kappa, n)

% alpha = circ_vmrnd(theta, kappa, n)

%   Simulates n random angles from a von Mises distribution, with preferred 

%   direction thetahat and concentration parameter kappa.

%

%   Input:

%     theta   - preferred direction

%     kappa   - width

%     n       - number of samples (default: 1)

%

%   Output:

%     alpha   - samples from von Mises distribution

%

%   Theta and kappa can be scalars, vectors, or matrices. If both or not

%   scalar, their dimensions should be identical and n should be 1. 

%

%   Examples

%    Y = circ_vmrnd(0,5)            - draw a sample from VM(0,5)

%    Y = circ_vmrnd(0,5,100)        - draw 100 samples from VM(0,5)

%    Y = circ_vmrnd(0,[5 10],100)   - draw 100 samples from VM(0,5) and 100 from VM(0,10)

%    Y = circ_vmrnd([-pi/2 pi/2],5) - draw one sample from VM(-pi/2,5) and one from VM(pi/2,5)

%

%   References:

%     Statistical analysis of circular data, Fisher, sec. 3.3.6, p. 49

%

% Circular Statistics Toolbox for Matlab



% By Philipp Berens and Marc J. Velasco, 2009

% Modified by Ronald van den Berg, 2011 



% input checking

warning off

if ~exist('n','var') 

    n=1;

end

if n==0

    Y = [];

    return

end



% handle all cases

if numel(kappa)==1 && numel(theta)==1     % both inputs are scalar

    input_dims = size(kappa);

    kappa=repmat(kappa,1,n);

    theta=repmat(theta,1,n);

elseif numel(kappa)==1 && numel(theta)>1  % kappa is scalar, theta is matrix

    input_dims = size(theta);

    theta = theta(:)';

    theta = repmat(theta,1,n);    

    kappa = ones(size(theta))*kappa;

elseif numel(kappa)>1 && numel(theta)==1  % kappa is matrix, theta is scalar

    input_dims = size(kappa);

    kappa = kappa(:)';

    kappa = repmat(kappa,1,n);    

    theta = ones(size(kappa))*theta;

elseif numel(kappa)>1 && numel(theta)>1   % both inputs are matrices

    if n>1

        error('Can only have n>1 when theta and kappa are scalars or vectors');        

    end

    if ~isequal(size(theta),size(kappa))

        error('Invalid input dimensions. If both theta and kappa is a vector or matrix, their dimensions should be the same');        

    end

    input_dims = size(kappa);

    kappa = kappa(:)';

    theta = theta(:)';

end



% use code from original circ_vmrnd to draw samples

a = 1 + sqrt((1+4*kappa.^2));

b = (a - sqrt(2*a))./(2*kappa);

r = (1 + b.^2)./(2*b);

valid = zeros(1,length(kappa));

z = zeros(size(kappa));

f = zeros(size(kappa));

c = zeros(size(kappa));

while ~all(valid)

    u(:,~valid) = rand(3,sum(~valid));    

    z(~valid) = cos(pi*u(1,~valid));

    f(~valid) = (1+r(~valid).*z(~valid))./(r(~valid)+z(~valid));

    c(~valid) = kappa(~valid).*(r(~valid)-f(~valid));       

    valid = u(2,:) < c .* (2-c) | ~(log(c)-log(u(2,:)) + 1 -c < 0);               

end

Y = theta + sign(u(3,:) - 0.5) .* acos(f);

Y = angle(exp(1i*Y));



% if kappa is very small, draw from a uniform distribution (the above method gives NaN for very small kappa's)

Y(kappa<1e-6) = 2*pi*rand(sum(kappa<1e-6),1)-pi;



% reshape back to original dimensions of input (add a dimension when n>1)

if n>1

    Y = reshape(Y,[max(input_dims) n]);

else

    Y = reshape(Y,input_dims);

end

warning on
