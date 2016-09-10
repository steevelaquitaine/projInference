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


%% 131004 Simulation of data to quantify the effect of prior strength on behavior.
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



%% 131005 Simulate fake data for 3 coherences and 4 priors with arbitrary widths.
%I choose values of k that produces simulations that look like the data.
%duration 15s for 50000 trials.
tic
k=[350 24 7.5 4.3 8.9 14.5 33];%true
[upo,d,coh,pstd]=simulateData(k,7000);
toc


%% 131008 Calculate logL of getting a data set given Girshick's model (26s).
tic
load data
kli(coh==0.24)=k(1);
kli(coh==0.12)=k(2);
kli(coh==0.06)=k(3);
kpi(pstd==80)=k(4);
kpi(pstd==40)=k(5);
kpi(pstd==20)=k(6);
kpi(pstd==10)=k(7);
[logL,ML]=GirshickML(upo,d,kpi,kli);
-sum(logL)
toc


%% 131008 Calculate logL of getting a data set given Girshick's model (26s)
%and looping over parameters.
load data
tic
kp1=2; %checked and it's ok
kp2=0:0.1:9.9;
kp3=0:2:18;
kp4=0:2:18;
kli(coh==0.24)=k(1);
kli(coh==0.12)=k(2);
kli(coh==0.06)=k(3);
kpi(pstd==80)=kp1;
tic
for j=1:1:100
    for k=1:10
        for l=1:10
            kpi(pstd==40)=kp2(j);
            kpi(pstd==20)=kp3(k);
            kpi(pstd==10)=kp4(l);
            [logL,ML]=GirshickML(upo,d,kpi,kli);
            slogL(i,j,k)=-sum(logL);
            
            hold all
            h=plot(1,1);
            plot(h,slogL(i,j,k),'ko')
            drawnow
        end
        toc
    end
end
figure;
[x,y,z]=meshgrid(1:10,1:10,1:10);
scatter3(x(:),y(:),z(:),100,slogL(:))


%% 131016 Look at the -sum(log likelihood) of parameters neighboring the true 
%parameters. The data are simulated data which we know the true parameters.
%We expect to see min(-sum(log likelihood)) lying at the true parameters.
%duration:1500 seconds
load data
getloglnearbyK

%% 131022 draw fit parameters & their std
load fit02
Sp=drawPrior(fitP.p,stdPa,[0.74559 2.77 8.74 33.25],[.24 .12 .06]);


%% 131023 Convolve llh with motor noise
%when motor noise is a von mises.
llh=vmPdfs(1:1:360,360,2.3);
Mn=vmPdfs(1:1:360,0,2.3);
cconv1=circConv(llh,Mn);
%when motor noise is a flat.
figure
llh=vmPdfs(1:1:360,360,2.3);
Mn=vmPdfs(1:1:360,0,0);
cconv2=circConv(llh,Mn);
%there is no motor noise.
figure
llh=vmPdfs(1:1:360,360,2.3);
Mn=vmPdfs(1:1:360,0,600);
cconv3=circConv(llh,Mn);


%% Prove that the maximum of the likelihood density is equivalent to the 
%probability of mi from the measurement density. 
%Note that P(mi|displayed) is the same in the measurement and llh
%densities and the llh is the measurement density shifted such that it now
%peaks at mi instead of the displayed direction.
%llh is the likelihood of displayed directions given measurement mi.
close all
mi=45;
mpdf=vmPdfs(1:1:360,180,3,'norm');
PmgivenDi=mpdf(mi)

llh=vmPdfs(1:1:360,mi,3,'norm');
PmigivenDi=llh(180)

hold all
plot(mpdf,'r')
plot(llh,'b')


%% 131021 find reasonable initial values for fitting by looking at mean data
clf
load data01
k=[250 24 7.5 4.3 8.9 14.5 33];
pred=drawDataFit(data,disp,coh,pstd,k,7000);



%% 131021 find reasonable initial values for fitting
%a)Try reasonable initial values based on qualitative fitting (eyes)
%(sub01)
clf
load data04
k=[250 34 7.5 8 14.5 28 50 1];
pred=drawDataFitwithFllh(data,disp,coh,pstd,k,7000);


%% 131022 find initial values for fitting
%with flat motor noise (or rather random estimates)
close all; clear
load data03
k=[18.1448982963345 3.45420350784409 0.0509269877940146 1.7018627128253 2.77114246357022 8.74031082097476 33.2498973961345 0.0325068906565146];%sub03-fit
drawRawDataFit(data,disp,coh,pstd,k,numel(data));

%% 131022 find initial values for fitting
%Full Bayesian inference model
initFullBI
simFullBI(data,disp,coh,pstd,k);

%% 131026 find reasonable initial values for fitting by looking at raw data
%Bayesian inference with changing likelihood model
initChangLLH
simChangLLH(data,disp,coh,pstd,k);

