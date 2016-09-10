%SLwjmPrefitInitpSub5to7onSuze.m
%
% author: steeve laquitaine
%   date: 150730 last modfi 150731
%purpose: fit WJM model to motion direction estimation data
%
%  usage:
%
%       o = SLwjmPrefitInitpSub5to7onSuze({'sub05','sub06','sub07'},[2 2 2 2.3 2.3 2.3 2.3 1 3],'experiment','vonMisesPrior') 
%    
%
%note: - tried to define prior on 1:1:360 but memory load on laptop. 
%      - in theory fit last 56 min per iteration with 5e3 trials and prior
%        defined for 100 values , no parfor.
%      - 4 minutes one iteration of the code (5e3 trials, prior defined for 100 directions, 
%        with parfor 10 cores (lapt).
%
%      - duration 2583.6sec (42 min,10 cores, 1 iter, 10 funEval, 5000 trials, 328 dots)
%      - duration 12651 sec (3 hours,10 cores, 15 iter, 150 funEval, 5000 trials, 328 dots)

function o = SLwjmPrefitInitpSub5to7onSuze(subjects,initp,varargin)

%delete(gcp)
%parpool('local')

%copy .m file 
%------------
mfilename = SLgetActivemFile;
myMfile = SLbackupMfileInWSpace(mfilename);
o.mfilename  = mfilename;
o.myMfile = myMfile;

%project path and data
%---------------------
fprintf('%s \n','(SLwjmPrefitInitpSub5to7onSuze) Gathering data ...')
fprintf('%s \n','(SLwjmPrefitInitpSub5to7onSuze) Please set dataPath ...')
dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');

if ~isempty(dataPath)
    cd(dataPath)
    fprintf('(SLwjmPrefitInitpSub5to7onSuze) Data path was set to : \n')
    fprintf(['(SLwjmPrefitInitpSub5to7onSuze) "',dataPath,'" \n'])
else
    fprintf('(SLwjmPrefitInitpSub5to7onSuze) Data path was not set...aborting \n')
    keyboard
end

varg = [varargin,'dataPath',dataPath];      %experiment and path
initp = repmat({initp},length(subjects),1); %init parameters

%load many subjects and create databases
parfor sub = 1 : length(subjects)
    
    fprintf('(SLwjmPrefitInitpSub5to7onSuze) Getting the subjects data...\n')
    
    %unique conditions
    db                                    = SLMakedatabank(subjects(sub),varg);
    db.estimatesDeg(db.estimatesDeg==0)   = 360;
    data{sub}                             = db.estimatesDeg;
    [taskC{sub},~,posC{sub}]              = SLuniqpair([db.Pstd  db.stimStrength  db.stimFeatureDeg]);
    
    %initialize model parameters for the 202 task conditions
    %[3 coh, 4 priors, random estimation and motor noise]
    %[0.24 0.12 0.06 80 40 20 10 prand km]
    pm0 = taskC{sub};
    pm0(taskC{sub}==0.24)   = initp{sub}(1);
    pm0(taskC{sub}==0.12)   = initp{sub}(2);
    pm0(taskC{sub}==0.06)   = initp{sub}(3);
    pm0(taskC{sub}==80)     = initp{sub}(4);
    pm0(taskC{sub}==40)     = initp{sub}(5);
    pm0(taskC{sub}==20)     = initp{sub}(6);
    pm0(taskC{sub}==10)     = initp{sub}(7);
    prand{sub}              = initp{sub}(8);
    km{sub}                 = initp{sub}(9);
    sCirad{sub}             = SLde2r(taskC{sub}(:,3),0); %s (known by subj,fixed) in rad
    cCi{sub}                = taskC{sub}(:,2);           %coh known by subj (fixed)
    kappaCi{sub}            = pm0(:,2);                  %sensory strength (free)
    kappa_sCi{sub}          = pm0(:,1);                  %prior strength (free)
    
    fprintf('%s %i %s \n','(SLwjmPrefitInitpSub5to7onSuze) - subjects ',sub,' ...done')
    
end

%search best parameters for many sub. in parallel
%-----------------------------------------------
tic
for sub = 1 : length(subjects); 
    savedfNames{sub} = ['dataCompleteSuzeSLwjmFitInitp1000iterTolX01sub0' num2str(sub) '.mat']; 
end
options = optimset('MaxIter',1000,'MaxFunEvals',10000,'TolFun',0.5,'TolX',0.1,'Display','iter');
fprintf('%s \n','count  fitp  logl')
count = 0;

%fit
parfor sub = 1 : length(subjects)
    
    fprintf('%s %i \n','(SLwjmPrefitInitpSub5to7onSuze) Fitting subject ',sub)
    fitPtmp         = initp{sub}; %initp
    datatmp         = data{sub};
    stmp            = sCirad{sub};
    ctmp            = cCi{sub};
    kappaCtmp       = kappaCi{sub};
    kappa_sCtmp     = kappa_sCi{sub};
    prandtmp        = prand{sub};
    kmtmp           = km{sub};
    posCtmp         = posC{sub};
    initptmp        = initp{sub};
    
    [fitPtmp,neglogltmp,exitflagtmp,outputfittmp] = fminsearch(@(fitPtmp) SLgetLogL_WJM(datatmp,...
        stmp,ctmp,kappaCtmp,kappa_sCtmp,prandtmp,kmtmp,fitPtmp,posCtmp,sub,count),...
        initptmp,options);
    
    fitP{sub}       = fitPtmp;          %best params
    neglogl(sub)    = neglogltmp;       % - logl
    exitflag{sub}   = exitflagtmp;      %info
    outputfit{sub}  = outputfittmp;
    
    slparsave(savedfNames{sub},o) %save

    fprintf('%s \n','(SLwjmPrefitInitpSub5to7onSuze) Fitting were successful...done')
    
end
tocs = toc;

%output
o.neglogl        = neglogl;
o.fitP           = fitP;
o.exitflag       = exitflag;
o.outputfit      = outputfit;
o.taskConditions = taskC;
o.fitDuration    = tocs;


%one fminsearch iteration
function [negloglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC,sub,count)

tic 

%count funEvals
count = count + 1;

%Simulating measurements x
Ndots   = 328;                            % n of dots
Ntrials = 5000;                           % n of trials

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
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2);              %proba
    PestGivModel(Ci,:) = (1-prand)*PestGivModel(Ci,:) + prand*(1/360);               %add random estimation
    PestGivModel(Ci,:) = SLcircConv(PestGivModel(Ci,:),vmPdfs(0:1:359,0,km,'norm')');%add motor noise

end

%get logl
for Ci = 1 : nCon

    PestGivModel(Ci,PestGivModel(Ci,:)<=0) = 10^-320;                   %no <0 due to convolution
    PestGivModel(Ci,:) = PestGivModel(Ci,:)./sum(PestGivModel(Ci,:),2); 
    logl(Ci) = sum(log(PestGivModel(Ci,data(posC==Ci))));
    
end

negloglAll = - sum(logl); %minimize -logl maximizes logl
fitP = fp;                %best fitp
tmpData = ['dataSuzeThisFunEvalSLwjmFitInitp1000iterTolX01sub0' num2str(sub)]; %save

o.negloglAll = negloglAll;
o.fitP       = fitP;
slparsave(tmpData,o)
tocs = toc;

fprintf('%i  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f %.2f  \n',count,fp,negloglAll,tocs)


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