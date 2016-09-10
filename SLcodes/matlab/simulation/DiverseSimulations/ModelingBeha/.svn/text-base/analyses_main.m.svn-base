%% Prior direction distribution biases perception as predicted by Bayesian inference
% INTRODUCTORY TEXT: 
    % Steeve laquitaine
    % started 12/09/07 - now
    % usage: this script summarizes all the analyses performed

%% figure 2. Experiment
%4 priors (10,20,40,80 degrees of std)
figure('color','w')
hold all
load('steeve_exp12_metho_Pstd080_mean225_coh006012024_dir36_t106_075_034perCoh_130217')
bar(5:10:355,task.parameter.dir.count/sum(task.parameter.dir.count)*100,'facecolor',[0.5 0 0],'edge','none')
load('steeve_exp12_metho_Pstd040_mean225_coh006012024_dir36_t101_075_031perCoh_130217')
bar(5:10:355,task.parameter.dir.count/sum(task.parameter.dir.count)*100,'facecolor',[1 0.2 0],'edge','none')
load('steeve_exp12_metho_Pstd020_mean225_coh006012024_dir36_t107_075_033perCoh_130217')
bar(5:10:355,task.parameter.dir.count/sum(task.parameter.dir.count)*100,'facecolor',[1 0.6 0],'edge','none')
load('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_073_033perCoh_130217')
bar(5:10:355,task.parameter.dir.count/sum(task.parameter.dir.count)*100,'facecolor',[0.75 0.75 0],'edge','none')
load('steeve_exp12_metho_Pstd010_mean225_coh006012024_dir36_t107_073_033perCoh_130217')
bar(5:10:355,task.parameter.dir.count/sum(task.parameter.dir.count)*100,'facecolor',[0.75 0.75 0],'edge','none')
xlabel('Displayed directions (degrees)')
ylabel('Probability (%)')


%% figure 3. *behavioral analyses*
tic
datlist.subjects={'sub01'};
fig.nm='steeve_descriptive_analyses';
[figOut,GenStat]=analyses(datlist,{'coh','Pstd','sample_dir'},fig);
toc


%% figure 4. *Modeling*
%30 min+
tic
datlist.subjects={'sub01'};
fig.nm='steeve_modeling';
% matlabpool 3
[data,disp,coh,pstd,Sp,pred,fitP,fitPbkp,R2,sdata,Fbp,F,stdPa,fitPt,logl,logL_bkp]=modelBehavior01(datlist,...
{'coh','Pstd','sample_dir'},... 
    fig,...
   [],[]);
toc

%% test without beta at all
%30 min+
tic
datlist.subjects={'sub01'};
fig.nm='steeve_modeling';
% matlabpool 10
[data,disp,coh,pstd,Sp,pred,fitP,fitPbkp,R2,sdata,Fbp,F,stdPa,fitPt,logl,logL_bkp]=modelBehavior01noBeta(datlist,...
{'coh','Pstd','sample_dir'},... 
    fig,...
   [],[]);
toc

%% model comparison
modelComparison

%% draw prior strength full Bayes
%131022 draw fit parameters & their std
load('Modeling/fitfullBayessub01.mat')
fitP01=fitP.p';
load('Modeling/fitfullBayessub02.mat')
fitP02=fitP.p';
load('Modeling/fitfullBayessub03.mat')
fitP03=fitP.p';
load('Modeling/fitfullBayessub04.mat')
fitP04=fitP.p';
fitP.p=[fitP01 fitP02 fitP03 fitP04];
s.nm={'s01','s02','s03','s04'};
s.num=size(fitP.p,2);
stdPa=zeros(size(fitP.p));%if stdPa was not calculated
truePrior=[0.74559 2.77 8.74 33.25];
coh=[.24 .12 .06];
Sp=drawPrior(fitP.p,stdPa,truePrior,coh,s);

%% draw prior strength changing lLH
%131022 draw fit parameters & their std
load('Modeling/fitchangLLHsub01.mat')
fitP01=fitP.p';
load('Modeling/fitchangLLHsub02.mat')
fitP02=fitP.p';
load('Modeling/fitchangLLHsub03.mat')
fitP03=fitP.p';
load('Modeling/fitchangLLHsub04.mat')
fitP04=fitP.p';
fitP.p=[fitP01 fitP02 fitP03 fitP04];
s.nm={'s01','s02','s03','s04'};
s.num=size(fitP.p,2);
stdPa=zeros(size(fitP.p));%if stdPa was not calculated
truePrior=[0.74559 2.77 8.74 33.25];
coh=[.24 .12 .06];
Sp=drawPrior(fitP.p,stdPa,truePrior,coh,s);

%% draw prior strength Heuristic
%131022 draw fit parameters & their std
load('fitHeuristicSub01_2.mat')
fitP01=fitP.p';
load('fitHeuristicSub02_2.mat')
fitP02=fitP.p';
load('fitHeuristicSub03_2.mat')
fitP03=fitP.p';
load('fitHeuristicSub04_2.mat')
fitP04=fitP.p';
fitP.p=[fitP01 fitP02 fitP03 fitP04];
s.nm={'s01','s02','s03','s04'};
s.num=size(fitP.p,2);
stdPa=zeros(size(fitP.p));%if stdPa was not calculated
truePrior=[0.74559 2.77 8.74 33.25];
coh=[.24 .12 .06];
Sp=drawPrior(fitP.p,stdPa,truePrior,coh,s);


