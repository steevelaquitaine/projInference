

%slfitRhoTauSigma.m
%
%
%
% author : steeve laquitaine
%purpose : fit the pp model by searching the best rho, tau and sigma
%          using maximum likelihood fitting with fminsearch Nelder-Mead



function [rho,tau,sigma,nglogl] = slfitRhoTauSigma(eta,nu)

%# of voxels
Nv = size(nu,1);

%I tested several scaling gain to the parameter to get better fit convergence 
%and speed up the fit and adding 20 (step size of 1 instead of the default 
%parameter magnitude *5%) produced lower log likelihood
%fp0 = [0.5; ones(Nv,1)*0.7; 0.3]+20;
% fp0 = [0.5; 0.05*randn(Nv,1)+0.7; 0.3]+20;
%update rho with stepsize 0.1 (+2 --> 5% of 2.1)
fp0(1) = 0.1 + 2*100;
%update tau with stepsize 0.1 (+1.4)
fp0(2:Nv+1,1) = ones(Nv,1)*0.7 + 1.4*100;
%update sigma with stepsize 0.1 (+1.8)
fp0(Nv+2) =  0.3 + 1.8*100;

%options were tuned by trial and error
options = optimset('maxIter',100000,'maxFunEvals',10000000,'TolFun',1,...
    'TolX',1e-1,'PlotFcns',@optimplotfval);

%to record the logl per iterations
global history; 
options.OutputFcn = slrecordfit;

%fit using maximum likelihood
[fp,nglogl,~,o] = fminsearch(@(fp) slgetloglppmodelStep2(eta,nu,fp),fp0,options);
%fp = fp - 20;
%get back unscaled values
fp(1) = fp(1) - 2*100;
fp0(2:Nv+1) = fp0(2:Nv+1) - 1.4*100;
fp0(Nv+2) =  fp0(Nv+2) - 1.8*100;


%plot logl convergence
figure('color','w'); plot(history.iter,history.objf)

%print convergence status
fprintf('%s \n',o.message)

%get best parameters
rho = fp(1);
tau = fp(2:end-1);
sigma = fp(end);

