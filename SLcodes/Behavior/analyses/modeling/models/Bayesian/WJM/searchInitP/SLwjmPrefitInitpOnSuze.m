%SLwjmPrefitInitpOnSuze.m
%
% author: steeve laquitaine
%   date: 150730 last modfi 150730
%purpose: fit WJM model to motion direction estimation data
%
%  usage:
%
%       o = SLwjmPrefitInitpOnSuze({'sub02'},[2 2 2 2.3 2.3 2.3 2.3 1 3],'experiment','vonMisesPrior') 
%    
%
%note: - tried to define prior on 1:1:360 but memory load on laptop. 
%      - in theory fit last 56 min per iteration with 5e3 trials and prioro = SLwjmPrefitInitpOnSuze({'sub02'},[2 2 2 2.3 2.3 2.3 2.3 1 3],'experiment','vonMisesPrior') 
%        defined for 100 values , no parfor.
%      - 4 minutes one iteration of the code (5e3 trials, prior defined for 100 directions, 
%        with parfor 10 cores (lapt).
%
%      - duration 2583.6sec (42 min,10 cores, 1 iter, 10 funEval, 5000 trials, 328 dots)
%      - duration 12651 sec (3 hours,10 cores, 15 iter, 150 funEval, 5000 trials, 328 dots)

function o = SLwjmPrefitInitpOnSuze(subjects,initp,varargin)

delete(gcp)
parpool('local')

%backup this .m file in the worspace in output variable
%------------------------------------------------------
mfilename = SLgetActivemFile;
myMfile = SLbackupMfileInWSpace(mfilename);
o.mfilename  = mfilename;
o.myMfile = myMfile;

%set path and data (ONCE)
%-------------
fprintf('%s \n','(SLfitBayesianModel_WJM) Gathering data ...')
fprintf('%s \n','(SLfitBayesianModel_WJM) Please set dataPath ...')
dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');
cd(dataPath)
varg = [varargin,'dataPath',dataPath]; %experiment  and path
db = SLMakedatabank(subjects,varg);

%get L(data|model) (ONCE)
%-----------------
data         = db.estimatesDeg;      
pstd         = db.Pstd;           
StimStrength = db.stimStrength;   
s            = db.stimFeatureDeg; 
[taskC,~,posC] = SLuniqpair([pstd StimStrength s]); %unique conditions

%initialize model parameters for the 202 task conditions (ONCE)
%3 coh, 4 priors, random estimation and motor noise
%[0.24 0.12 0.06 80 40 20 10 prand km]              
pm0 = taskC;
pm0(taskC==0.24)   = initp(1);
pm0(taskC==0.12)   = initp(2);
pm0(taskC==0.06)   = initp(3);
pm0(taskC==80)     = initp(4);
pm0(taskC==40)     = initp(5);
pm0(taskC==20)     = initp(6);
pm0(taskC==10)     = initp(7);
prand              = initp(8);
km                 = initp(9);
sCirad    = SLde2r(taskC(:,3),0); %s (known by subj,fixed) in rad
cCi       = taskC(:,2);           %coh known by subj (fixed)
kappaCi   = pm0(:,2);             %sensory strength (free)
kappa_sCi = pm0(:,1);             %prior strength (free)

%search best parameters
%----------------------
tic
options = optimset('MaxIter',1000,'MaxFunEvals',10000,'TolFun',0.5,'TolX',0.1,'Display','iter');
[fitP,neglogl,exitflag,outputfit] = fminsearch(@(fitP) SLgetLogL_WJM(data,sCirad,cCi,kappaCi,kappa_sCi,prand,km,fitP,posC),...
    initp,options);
tocs = toc;

%output
o.neglogl = neglogl;
o.fitP = fitP;
o.exitflag = exitflag;
o.outputfit = outputfit;
o.taskConditions = taskC;
o.fitDuration = tocs;
save('dataSLwjmPrefitInitp1000iterTolX01')


%one fminsearch iteration
function [negloglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC)

%Simulating measurements x
Ndots = 328;                              % number of dots
Ntrials = 5000;                           % number of trials

%prior parameters
Ns = 100;
svec = linspace(1,360,Ns+1)';
svec = SLde2r(svec(1:end-1),0);  % vector of hypothesized motion directions

%predicted estimates space
predEst = 1:1:360;
nCon = max(posC);

%scale prand (for tolX, speed up)
prand = prand./1e6;
km = km*1e30;

parfor Ci = 1 : nCon
        
    %measurements
    coherent = (rand(Ntrials,Ndots) < c(Ci)); % which dots are moving coherently
    x = circ_vmrnd(s(Ci)*ones(Ntrials,Ndots),kappa(Ci));
    x = coherent .* x + (1-coherent) .* rand(Ntrials,Ndots)*2*pi; % Ntrials by Ndots
    
    %prior
    Mean_s = SLde2r(225,0);          % assumed known to subj (fixed)
    prior = exp(kappa_s(Ci) * cos(svec - Mean_s));
    
    % Likelihoods L(s) = p(x|s)
    x = permute(repmat(x, [1 1 Ns]),[3 1 2]); % to make dimensions match: Ns by Ntrials by Ndots
    like_per_dot = (1-c(Ci))/2/pi + c(Ci)/2/pi/besseli(0,kappa(Ci)) * exp(kappa(Ci) * cos(bsxfun(@minus,x,svec))); %(1-c)/2pi + c*V(svec,x,kappa)
    like = prod(like_per_dot,3);
    
    % Posteriors p(s|x)
    posterior = bsxfun(@times, prior, like); % Ns by Ntrials
    
    % Estimation: posterior means (circular)
    xtemp = sum(bsxfun(@times, cos(svec), posterior));
    ytemp = sum(bsxfun(@times, sin(svec), posterior));
    posteriormeans = atan2(ytemp,xtemp); % 1 by Ntrials)
    posteriormeans = SLra2d(posteriormeans);
    
    %estimate distributions
    PestGivModel(Ci,:) = hist(posteriormeans,predEst);
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    PestGivModel(Ci,:) = (1-prand)*PestGivModel(Ci,:) + prand*(1/360);               %add random estimation
    PestGivModel(Ci,:) = SLcircConv(PestGivModel(Ci,:),vmPdfs(0:1:359,0,km,'norm')');%add motor noise

    %fprintf('%s %i \n','Cond: ',Ci)
end

%get logl
for Ci = 1 : nCon

    PestGivModel(Ci,PestGivModel(Ci,:)<=0) = 10^-320;                   %get rid of negatives due to convolution
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); %proba
    logl(Ci) = sum(log(PestGivModel(Ci,data(posC==Ci))));
    
    %fprintf('%s %.0f \n','logl: ',logl)

end
negloglAll = - sum(logl); %minimize -logl maximizes logl
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