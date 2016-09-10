clc
%List of simulations & modeling analyses


%% Girshick estimator.
%Simulate estimates' densities for a couple of km & kp when we are
%away or close to prior's mean. Estimate densities' mean & std.
k.m=1; 
k.p=0; 
[mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs,means,stds]=GirshickEstimator(k);



%% Range of estimate predictions of Girshick estimator.
%To see the estimates' densities predicted byy girshick estiator model.
[mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs]=GirshickPredictions;



%% 131004 Simulation to calculate standard deviation of von Mises and correspondence
%to gaussians
%Data derived estimates of the k parameters are [25 10 2.5 0.74559 2.77 8.74 33.25];
%see what the distributions look like (with polar too)
%The best gaussian approximation of a von mises is when the mean of the 
%density=180.
figure('color','w')
[S1,S2]=simulateVMStd(1:1:360,180,100);
S2


%% 131004 simulation of data to quantify the effect of prior strength on behavior.
%km=2 (obtained from model fitting). 
%Simulated estimates don't change anymore when prior is strong.
%note: Larger noise when prior is weak makes it difficult 
%to quantify the changes.
figure('color','w'); hold all
title('Predictions do not change when prior is strong','fontsize',14)
c=colormap('lines');

klset=[0 16.5*(1:1:5)];
for i=1:numel(klset); 
    step(i)=klset(i);%i*5-5; 
    [d(:,i),o(:,i)]=simulateBasicData(33,step(i),c(i,:)); 
end
drawPublishAxis

%quantify change
figure('color','w'); hold all
deltao=mean(abs(o-225));
plot(step(1:end),deltao,'-o')
ylabel('Absolute change in simulated estimates:mean(abs(es-225));degrees)')
xlabel('Prior strength (k_p)')
drawPublishAxis



%% 131005 simulate data for 3(/4) coherences(/priors)
initsimData
[upo,d,coh,pstd]=simulateData(k,7000);


%% 131008 Calculate logL of data given Girshick's model
%get logl given a set of fit parameters
getGirshickLogL


%% 131008 calculate logL of data given Girshick's model
%get logl from a grid of fit parameters
getGirshickLogLgrid


%% 131016 Look at the -sum(log likelihood) of parameters neighboring the true 
%parameters. The data are simulated data which we know the true parameters.
%We expect to see min(-sum(log likelihood)) lying at the true parameters.
%duration:1500 seconds
load data
getloglnearbyK


%% 131023 draw llh convolved with different motor noises
seeConvMotorNoises


%% Prove that P(mean of llh) = P(mi)
%P(mi|displayed) is the same in the measurement and llh densities
%and llh is measurement density shifted such that it now peaks at 
%mi instead of the displayed direction. llh is the likelihood of displayed 
%directions given measurement mi.
ProvePllhmeanPmi


%% 131022. find initial values for fitting
%Full Bayesian inference model
initFullBI
simFullBI(data,disp,coh,pstd,k);

%% 131123. draw raw data and simulate data from model of "CHANGING LIKELIHOOD"
%We use it also to model full Bayes inference
%to find best initial parameters
tic
%parpool 12
initChangLLH
[pred,Kinput,Kreal]=simChangLLH(data,disp,coh,pstd,k,'TheoreticalMean',200);
toc

%% 131202. check fitting procedure for "CHANGING LIKELIHOOD" model
checkfitchangLLH


%% 131214. draw simulations of the softmax competition stage of the "heuristic" model (16s)
initHeuristicSoftmax
[Ppriormn,Pllhmn,Prnew]=simHeuristicSoftmax(kp,kml,Pr,Beta);


%% 131212. draw predictions of "heuristic" model (16s)
%to find best initial parameters
tic
%parpool 12
initHeuristic
[pred,Kinput,Kreal,mixtureProba]=simHeuristic(data,disp,coh,pstd,k,'TheoreticalMean',1);
toc


%% 131216. check fitting procedure for "Bayes-like heuristic" model
checkFitHeuristic


%% 131212. draw predictions of "extended" Bayesian model (16s)
%to find best initial parameters
tic
%parpool 12
%close all
initExtendedBI
[pred,Kinput,Kreal]=simExtendedBI(data,disp,coh,pstd,k,'TheoreticalMean',1);
toc


%% 1312124. Compare a von mises density with a Gaussian density in case 
%of our priors
vonMisesVsGaussian



