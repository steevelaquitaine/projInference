%slFitWJM.m
%
%
% author: steeve laquitaine
%   date: 150821 updated 150826...currently writing
%purpose: fit WJM model to motion direction estimation data
%
%  usage:
%
%       o = slFitWJM({'sub01'},[2 2 2 2.3 2.3 2.3 2.3 1e3 3e-30],'experiment','vonMisesPrior') 
%       o = slFitWJM({'sub01'},[2 2 2 2.3 2.3 2.3 2.3 1 3],'experiment','vonMisesPrior') 
%
%
%-----
%note: 
%-----
%param setting
%   
%      - "prand" set at 1e3 (-> 1e-3 same as initp for switching)    
%      - "km" set at 3e-30  (-> 3 same as switching)
%
%Test
%      - the code has been tested with 10 dots, 100 trials (quick iter) and fminsearch
%        minimizes the -logl has expected.
%
%      - the current params ("prand" and "km") are same as for switching model
%        and also produce smooth predicted estimate densities.
%
%      - "prand" affects the obj.fun a lot (0 (1.6e+6) vs 1(87e+3)).
%        because 0 means some data (e.g., mistakes) are impossible based on the model. We
%        need non-zero prand (lapse rate).
%
%      - "km" affects the obj.fun a lot (0(31e3) vs 3(85e3)). This is
%        because convolution smoothes up model predicted estimates
%        densities which otherwise are very noisy (random simulation) 
%        and raises the obj. fun values 
%
%      - # of trials also reduces obj.fun a lot (100 (85e3) vs. 5000(74e3))
%        because they smooth up predicted estimate densities.
%
%   
%Duration
%      - optimization returns everytime <0 params is observed. 
%      - by setting TolX to 0.1, I speed up optimization termination.
%      - tried to define prior on 1:1:360 but memory load (laptop) so defined on 100 values. 
%      - in theory fit last 56 min per iter (5e3 trials, prior defined for 100 values, no parfor).
%      - 4 minutes per iteration of the code (5e3 trials, prior defined for 100 directions, 
%        with parfor 10 cores (lapt).
%
%      - duration 2583.6sec (42 min,10 cores, 1 iter, 10 funEval, 5000 trials, 328 dots)
%      - duration 12651 sec (3 hours,10 cores, 15 iter, 150 funEval, 5000 trials, 328 dots)
%
%      - running on Sherlock with parfor on 16 cores (10 cores, 1000 iter -> 8 days, 16 cores should be <7 days).
%        we choose 1000 iterations because it should last 7 days max (limit on Sherlock and Barley)
%
function o = slFitWJM(subjects,initp,varargin)

    %add code libraries
    addpath(genpath('~/proj/'))

    %delete(gcp)
    %local = parpool(16)
    %parpool('local')

    %----------------
    %backup this code 
    %----------------
    %mfilename    = SLgetActivemFile;
    %myMfile      = SLbackupMfileInWSpace(mfilename);
    %o.mfilename  = mfilename;
    %o.myMfile    = myMfile;

    %----------------------------------
    %get project path and organize data
    %----------------------------------
    fprintf('%s \n','(slFitWJM) Gathering data ...')
    fprintf('%s \n','(slFitWJM) Please set dataPath ...')
    dataPath = uigetdir(cd,'Pick a project e.g., /dataPsychophy/Exp01...');
    %dataPath = '~/data';

    if ~isempty(dataPath)
        cd(dataPath)
        fprintf('(slFitWJM) Data path set to : \n')
        fprintf(['(slFitWJM) "',dataPath,'" \n'])
    else
        fprintf('(slFitWJM) Data path not set...aborting \n')
        keyboard
    end

    varg = [varargin,'dataPath',dataPath];      %experiment and path
    initp = repmat({initp},length(subjects),1); %init parameters

    %load subject and create database
    sub = 1;
    fprintf('(slFitWJM) Getting subject data...\n')

    %conditions
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

    fprintf('%s %s %s \n','(slFitWJM) - Database for subjects ',subjects{sub},' ...organized')

    %----------------------    
    %fit WJM model to data
    %----------------------
    tic

    %saved filename
    for sub = 1 : length(subjects); 
        savedfNames{sub} = ['fitBackupFinal.mat']; 
    end

    %update
    fprintf('%s \n','-----------')
    fprintf('%s \n',' fitp  logl')
    fprintf('%s \n','-----------')

    %fit
    fprintf('%s %i \n','(slFitWJM) Fitting subject ',subjects{sub})
    %options = optimset('OutputFcn', @myoutput,'MaxIter',1000,'MaxFunEvals',10000,'TolFun',0.5,'TolX',0.1,'Display','iter');
    options = optimset('OutputFcn', @myoutput,'MaxIter',100,'MaxFunEvals',100,'TolFun',0.5,'TolX',0.1,'Display','iter');

    %init stored path and fit variables
    mkdir(['fitBackup' subjects{sub}])
    cd(['fitBackup' subjects{sub}])
    history.params    = [];
    history.iteration = [];
    history.funccount = [];
    history.fval      = [];

    %parameters
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

    %fit
    [fitPtmp,neglogltmp,exitflagtmp,outputfittmp] = fminsearch(@(fitPtmp) SLgetLogL_WJM(datatmp,...
        stmp,ctmp,kappaCtmp,kappa_sCtmp,prandtmp,kmtmp,fitPtmp,posCtmp,subjects{sub}),...
        initptmp,options);

    %backup
    fitP{sub}       = fitPtmp;          %best params
    neglogl(sub)    = neglogltmp;       % - logl
    exitflag{sub}   = exitflagtmp;      %info
    outputfit{sub}  = outputfittmp;
    slparsave(savedfNames{sub},o)       %save

    fprintf('%s \n','(slFitWJM) Fitting was successful...done')
    tocs = toc;

    %output
    o.neglogl        = neglogl;
    o.fitP           = fitP;
    o.exitflag       = exitflag;
    o.outputfit      = outputfit;
    o.taskConditions = taskC;
    o.fitDuration    = tocs;


    %fit output
    function stop = myoutput(x,optimvalues,state);

        stop = false;
        if isequal(state,'iter')
          history.params    = [history.params    ; x                    ];
          history.iteration = [history.iteration ; optimvalues.iteration];
          history.funccount = [history.funccount ; optimvalues.funccount];
          history.fval      = [history.fval      ; optimvalues.fval     ];
        end

        %save filename fname
        save('fitBackup','history')
    end

    %one fminsearch iteration
    function [negloglAll,fitP,PestGivModel] = SLgetLogL_WJM(data,s,c,kappa,kappa_s,prand,km,fp,posC,sub)

        tic 

        %sanity check
        if any(fp<0)
            
            negloglAll = +inf;
            fitP = [];
            PestGivModel = [];
            
            fprintf('(SLgetLogL_WJM) WARNING ! Unreasonable (<0) param value.')
            return
        end

        %Simulating measurements x
        Ndots   = 10;%328; %# dots
        Ntrials = 10;%5000; %# trials

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
            x = slCirc_vmrnd(s(Ci)*ones(Ntrials,Ndots),kappa(Ci));
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

        %backup
        o.negloglAll = negloglAll;
        o.fitP       = fitP;
        tocs = toc;

        fprintf('%.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f  %.2f %.2f  \n',fp,negloglAll,tocs)
    end    
end