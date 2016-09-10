
% slSimulateBayesianModel.m
%
%     author: steeve laquitaine
%       date: 140601 updated 150110
%    purpose: simulate motion direction estimation behavior with a Bayesian
%             observer (learnt priors, cardinal priors, motor noise, random estimation).
%
%
%inputs:
%
%       'dbBayesModel.mat' : is in the current directory and contains the
%                            experiment parameters (from SLmakeDatabank.m)
%
%options:
%
%   'MAPReadout'                  : estimate is the argmax of the posterior
%   'SamplingReadout'             : estimate is sampled from the posterior
%   'modelTheorecticalPredictions': simulation of calculated estimate mean and std.
%
%references:
%
%     -Hurliman et al, 2002,VR
%     -Stocker&Simoncelli,2006,NN
%     -Girshick&Simoncelli,2011,NN
%     -Chalk&Series,2012,JoV
%
%
%Usage:
%
%pathLab = '/Users/steeve/Steeve_psychophy_data/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior';
%pathHome = '/Users/steeve_laquitaine/Dropbox/Project_withJustin/data/psychophy_dataForPlayingAround/Exp01_Variance_Pstd010020040080_coh006012024_dirfull_Pmean225_fix10_mot035_fB01_speed28_PwMate_reBalanced_stat4_v2/MostUpdated/Code_004mainfinalCardinalsOutRndInitPosSymPrior'
% 
% slSimulateBayesianModel([10 8 4 1.74 4.77 10.74 34.25 0 0.001 15],...
%     'experiment','vonMisesPrior','MAPReadout');
%
% [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=slSimulateBayesianModel({'sub03'},...
%     [80 40 3 30 30 30 30 0 0.001 15],...
%     'dataPath',pathLab,...
%     'experiment','bimodalPrior');
%
% [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL]=slSimulateBayesianModel({'sub01'},...
%     [26 5 1.85 0.7 2.7 8.7 33 NaN 0 3000],...
%     'dataPath',pathLab,...
%     'experiment','vonMisesPrior',...
%     'MAPReadout',...
%     'modelTheoreticalPredictions');
% 
% slSimulateBayesianModel({'sub01'},...
%     [26 5 1.85 0.7 2.7 8.7 33 NaN 0 3000],...
%     'experiment','vonMisesPrior',...
%     'MAPReadout',...
%     'modelTheoreticalPredictions');

function [pred,coh,pstd,d,priorModes,simP,negElogL,negSumlogL] = slSimulateBayesianModel(simP,varargin)

%(case simulated estimate mean and std based on sampling
%-------------------------------------------------------

%load motion direction experiment parameters
ws         = load('dbBayesModel');
databank   = ws.databank;
d          = cell2mat(databank.data(:,(strcmp(databank.nm,'FeatureSample'))==1));
coh        = cell2mat(databank.data(:,(strcmp(databank.nm,'StimStrength'))==1));
pstd       = cell2mat(databank.data(:,(strcmp(databank.nm,'Pstd'))==1));
priorModes = cell2mat(databank.data(:,(strcmp(databank.nm,'priormodes'))==1));

%get fit parameters. They should be in this order.
%'coh24','coh12','coh6','pstd80','pstd40','pstd20','pstd10','kc','Prand','km'
%likelihood strengths
kl1 = simP(1); kl2 = simP(2); kl3 = simP(3);

%prior strengths
klearnt1 = simP(4); klearnt2 = simP(5); klearnt3 = simP(6); klearnt4 = simP(7); kc = simP(8);

%random error and motor noise
Prandom = simP(9);
km = simP(10);
ulall=1:1:360;

TheModel = 'withoutCardinal';
if isnan(kc); TheModel ='withCardinal';end
    

%Penalize parameter values out of range. Immediately go to the next
%iteration. It is important to have this here instead of having it at the
%end of the code to save processing time.
%Constrain other fraction parameters between 0 and 1, and other parameters
%as>0.
%and case a parameter that is not cardinal prior strength is missing (NaN)
if Prandom>1
    negElogL=inf;
    return
end
if any(simP<0)
    negElogL=inf;
    return
end
if any(isnan(simP([1,2,3,4,5,6,7,9,10])))
    negElogL=inf;
    fprintf('%s \n',['(getLogl) One of your fit parameter',...
        'that is not Kcardinal is NaN'])
end

%get the modes of the prior
SetOfpriorModes = SLuniqpair(priorModes);

%(case learnt prior is a von Mises)
%----------------------------------
%if the prior is a von Mises, the unique mode of the prior(s) are those
%found in the databank from (SLMakedatabank).
%In our unimodal experiment it should be [225].
if sum(strcmp(varargin,'vonMisesPrior'))==1
    priorShape='vonMisesPrior';
    
    %get priors modes (225)
    modesPrior1=SetOfpriorModes(1,:);
    modesPrior2=SetOfpriorModes(1,:);
    modesPrior3=SetOfpriorModes(1,:);
    modesPrior4=SetOfpriorModes(1,:);
    
    %find trials for each prior
    Prior1=pstd==80;
    Prior2=pstd==40;
    Prior3=pstd==20;
    Prior4=pstd==10;
end

%(case bimodal prior)
%--------------------
%get the four prior conditions in this order:
if sum(strcmp(varargin,'bimodalPrior'))
    priorShape='bimodalPrior';
    
    %get priors modes
    %[145 305]
    %165 285]
    %[185 265]
    %[205 245]
    modesPrior1=SetOfpriorModes(1,:);
    modesPrior2=SetOfpriorModes(2,:);
    modesPrior3=SetOfpriorModes(3,:);
    modesPrior4=SetOfpriorModes(4,:);
    
    %find trials for each prior
    Prior1=priorModes(:,1)==modesPrior1(1,1) & priorModes(:,2)==modesPrior1(1,2);
    Prior2=priorModes(:,1)==modesPrior2(1,1) & priorModes(:,2)==modesPrior2(1,2);
    Prior3=priorModes(:,1)==modesPrior3(1,1) & priorModes(:,2)==modesPrior3(1,2);
    Prior4=priorModes(:,1)==modesPrior4(1,1) & priorModes(:,2)==modesPrior4(1,2);
end

%-------
%SET LLH
%-------
%mean
%di=unique(d);
motdir=unique(d);

%strength
coh24=coh==0.24;
coh12=coh==0.12;
coh06=coh==0.06;


%1.Bayesian inference
%--------------------
%get a matrix (360 possible MAPs, 360 possible motion directions) of
%likelihood values for each data (MAPs,rows) and each motion directions
%(columns). We get this matrix for each of the 12 condition (3 coh x 4
%learnt priors) of the experiment. This is independent of subjects' data.
%I have checked what those matrices look like and the results are
%intuitifs.

%(case MAP readout)
%-----------------------
if sum(strcmp(varargin,'MAPReadout'))==1
    
    %get likelihood of data in condition with coherence 24
    [~,likeMAP11]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior1,klearnt1,kc,0,priorShape,TheModel,varargin);
    [~,likeMAP12]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior2,klearnt2,kc,0,priorShape,TheModel,varargin);
    [~,likeMAP13]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior3,klearnt3,kc,0,priorShape,TheModel,varargin);
    [~,likeMAP14]=SLGirshickBayesLookupTable(ulall,motdir,kl1,modesPrior4,klearnt4,kc,0,priorShape,TheModel,varargin);
        
    %get likelihood of data in condition with coherence 12
    [~,likeMAP21]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior1,klearnt1,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP22]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior2,klearnt2,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP23]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior3,klearnt3,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP24]=SLGirshickBayesLookupTable(ulall,motdir,kl2,modesPrior4,klearnt4,0,kc,priorShape,TheModel,varargin);
    
    %get likelihood of data in condition with coherence 6
    [~,likeMAP31]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior1,klearnt1,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP32]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior2,klearnt2,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP33]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior3,klearnt3,0,kc,priorShape,TheModel,varargin);
    [MAP,likeMAP34]=SLGirshickBayesLookupTable(ulall,motdir,kl3,modesPrior4,klearnt4,0,kc,priorShape,TheModel,varargin);
end

%(case Sampling readout)
%-----------------------
if sum(strcmp(varargin,'SamplingReadout'))==1
    
    %get likelihood of data in condition with coherence 24
    [~,likeMAP11]=SLBayesSamplingLookupTable(ulall,motdir,kl1,modesPrior1,klearnt1,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP12]=SLBayesSamplingLookupTable(ulall,motdir,kl1,modesPrior2,klearnt2,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP13]=SLBayesSamplingLookupTable(ulall,motdir,kl1,modesPrior3,klearnt3,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP14]=SLBayesSamplingLookupTable(ulall,motdir,kl1,modesPrior4,klearnt4,0,kc,priorShape,TheModel,varargin);
    
    %get likelihood of data in condition with coherence 12
    [~,likeMAP21]=SLBayesSamplingLookupTable(ulall,motdir,kl2,modesPrior1,klearnt1,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP22]=SLBayesSamplingLookupTable(ulall,motdir,kl2,modesPrior2,klearnt2,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP23]=SLBayesSamplingLookupTable(ulall,motdir,kl2,modesPrior3,klearnt3,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP24]=SLBayesSamplingLookupTable(ulall,motdir,kl2,modesPrior4,klearnt4,0,kc,priorShape,TheModel,varargin);
    
    %get likelihood of data in condition with coherence 6
    [~,likeMAP31]=SLBayesSamplingLookupTable(ulall,motdir,kl3,modesPrior1,klearnt1,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP32]=SLBayesSamplingLookupTable(ulall,motdir,kl3,modesPrior2,klearnt2,0,kc,priorShape,TheModel,varargin);
    [~,likeMAP33]=SLBayesSamplingLookupTable(ulall,motdir,kl3,modesPrior3,klearnt3,0,kc,priorShape,TheModel,varargin);
    [MAP,likeMAP34]=SLBayesSamplingLookupTable(ulall,motdir,kl3,modesPrior4,klearnt4,0,kc,priorShape,TheModel,varargin);
end

%Probability of getting a given estimate given Bayesian inference
%now get matrix 'PupoGivenBI' of likelihood values (upos=1:1:360,trials)
%for possible values of upo (rows) for each trial (column)
PupoGivenBI=nan(numel(MAP),numel(d));
for i=1:numel(motdir)
    
    %get displayed motion direction for this trial
    thisd=d==motdir(i);
    
    %get likelihood of data in condition with learnt prior 80
    PupoGivenBI(:,thisd&coh24&Prior1)=likeMAP11(:,motdir(i(ones(sum(thisd&coh24&Prior1),1))));
    PupoGivenBI(:,thisd&coh12&Prior1)=likeMAP21(:,motdir(i(ones(sum(thisd&coh12&Prior1),1))));
    PupoGivenBI(:,thisd&coh06&Prior1)=likeMAP31(:,motdir(i(ones(sum(thisd&coh06&Prior1),1))));
    
    %get likelihood of data in condition with learnt prior 40
    PupoGivenBI(:,thisd&coh24&Prior2)=likeMAP12(:,motdir(i(ones(sum(thisd&coh24&Prior2),1))));
    PupoGivenBI(:,thisd&coh12&Prior2)=likeMAP22(:,motdir(i(ones(sum(thisd&coh12&Prior2),1))));
    PupoGivenBI(:,thisd&coh06&Prior2)=likeMAP32(:,motdir(i(ones(sum(thisd&coh06&Prior2),1))));
    
    %get likelihood of data in condition with learnt prior 20
    PupoGivenBI(:,thisd&coh24&Prior3)=likeMAP13(:,motdir(i(ones(sum(thisd&coh24&Prior3),1))));
    PupoGivenBI(:,thisd&coh12&Prior3)=likeMAP23(:,motdir(i(ones(sum(thisd&coh12&Prior3),1))));
    PupoGivenBI(:,thisd&coh06&Prior3)=likeMAP33(:,motdir(i(ones(sum(thisd&coh06&Prior3),1))));
    
    %get likelihood of data in condition with learnt prior 10
    PupoGivenBI(:,thisd&coh24&Prior4)=likeMAP14(:,motdir(i(ones(sum(thisd&coh24&Prior4),1))));
    PupoGivenBI(:,thisd&coh12&Prior4)=likeMAP24(:,motdir(i(ones(sum(thisd&coh12&Prior4),1))));
    PupoGivenBI(:,thisd&coh06&Prior4)=likeMAP34(:,motdir(i(ones(sum(thisd&coh06&Prior4),1))));
end

%scale to probability.
Z_=sum(PupoGivenBI);
Z=Z_(ones(size(PupoGivenBI,1),1),:);
PupoGivenBI=PupoGivenBI./Z;

%get probabilities of the range of percepts "upo" given random estimation
%------------------------------------------------------------------------
PupoGivenRand=ones(360,numel(d))/360;

%calculate probability of percepts "upo" given the model
PBI=1-Prandom;
PupoGivenModel=PupoGivenBI*PBI + PupoGivenRand*Prandom;

%check PupoGivenModel sum to 1
if ~unique(sum(PupoGivenModel))==1
    fprintf('%s \n','Something s wrong. PupoGivenModel are probabilties and should sum to 1')
    keyboard
end

%convolve with motor noise
%------------------------------------------------------------------------
%Now we shortly replace upo=1:1:360 by upo=0:1:359 because motor noise
%distribution need to peak at 0 and vmPdfs function needs 'x' to contain
%the mean '0' to work. Then we set back upo to its initial value. This have
%no effect on the calculations.
%upo=1:1:360;
upo=0:1:359;
Pmot=vmPdfs(upo,0,km,'norm');
Pmotforconv=Pmot(:,ones(1,numel(d)));
PestimateGivenModel=SLcircConv(PupoGivenModel,Pmotforconv);

%get read of extremely small negative numbers that occurs because of
%instability numerical instability. Set them to zero.
PestimateGivenModel(PestimateGivenModel<0)=0;


%(case model's theortical predictions about estimates mean and std)
%------------------------------------------------------------------
if sum(strcmp(varargin,'modelTheoreticalPredictions'))==1
    
    %status
    fprintf('%s \n','(SLsimulateBayesianModel) Calculating model',...
        ' theortical predictions about estimates mean and std at convergence')
    
    %get predictions densities sorted per condition
    %We extract the predicted densities of estimates for each single
    %condition "cond" of the experiment.
    %cond(idxCond(thistrial),:) tells the condition (coh pstd displ) for each
    %trial SLuniqpair by default gives the index of the first trial with the
    %condition in "cond".
    %"PestimateGivenModelUniq"'s rows are possible estimates and columns are
    %experimental condition in rows of "cond".
    [cond,idxCondUniqtrial,~]=SLuniqpair([pstd coh d]);
    PestimateGivenModelUniq=PestimateGivenModel(:,idxCondUniqtrial);
    
    %sort everything to match the data sorting (ascending order)
    [cond,PosSorted]=sortrows(cond,[-2 -1]);
    PestimateGivenModelUniq=PestimateGivenModelUniq(:,PosSorted);
    
    %calculate circular mean and std or predictions for each cond (columns)
    %Use circular statistics
    meanPred=nan(size(PestimateGivenModelUniq,2),1);
    stdPred=nan(size(PestimateGivenModelUniq,2),1);
    for i = 1 : size(PestimateGivenModelUniq,2)
        data = SLcircWeightedMeanStd(MAP, PestimateGivenModelUniq(:,i));
        meanPred(i)= data.deg.mean;
        stdPred(i) = data.deg.std;
    end
    
    %draw models theoretical predictions
    drawMeanPreCentered(meanPred,stdPred,cond)
    
else
    
    %(case model's  predictions about estimates mean and std by averaging
    %trial samples)
    %---------------------------------------------------------------------
    
    %set upo to initial values any case we use it later
    upo=1:1:360;
    
    %normalize to probability; If we don't normalize, more random
    %choice always increase probability of observing data causing larger
    %probability of random choice to prevail. To avoid that we need to normalize
    %to probabilities that sum to 1. It also makes intuitive sense to deal
    %with probabilities.
    Z_=sum(PestimateGivenModel);
    Z=Z_(ones(size(PestimateGivenModel,1),1),:);
    PestimateGivenModel=PestimateGivenModel./Z;
    
    %trial-predictions (sample estimate density). Variability in estimate
    %should reflect variability in measurement density.
    for i=1:length(d)
        pred(i)=randsample(1:1:360,1,'true',PestimateGivenModel(:,i));
    end
    
    %draw sampled predictions
    %------------------------
    %draw predictions and data
    %case von Mises Prior
    %--------------------
    if sum(strcmp(varargin,'vonMisesPrior'))
        SLdrawCircStat(pred,d,coh,pstd);
        SLdrawCircStat(data,d,coh,pstd);
    end
    
    %case bimodal Prior
    %--------------------
    if sum(strcmp(varargin,'bimodalPrior'))
        SLdrawCircStat(pred,d,coh,priorModes(:,2)-priorModes(:,1));
        SLdrawCircStat(data,d,coh,priorModes(:,2)-priorModes(:,1));
    end
    
    %get loglikelihood of data
    %-------------------------
    %single trial's measurement, its position(row) for each trial(col) and its
    %probability (also maxlikelihood of trial's data). Checked many times. It
    %works.
    %make sure sub2ind inputs are the same size
    if sum(size(data)~=size((1:1:numel(d))'))==2
        data=data';
    end
    idx=sub2ind(size(PestimateGivenModel),data,(1:1:numel(d))');
    PdataGivenModel=PestimateGivenModel(idx);
    
    %We use log likelihood because likelihood is so small that matlab cannot
    %encode it properly (numerical unstability). We can use single trials log
    %likelihood to calculate AIC in the conditions that maximize differences in
    %predictions of two models.
    logL_pertrial=log(PdataGivenModel);
    
    %We use -sum(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    negSumlogL=-sum(log(PdataGivenModel));
    
    %We use E(log likelihood) as an objective function to minimize with
    %matlab non linear optimization search.
    negElogL=-nanmean(log(PdataGivenModel));
    
    %Look at fitting. It is 3X faster without drawing.
    fprintf('%.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.2f   %.05f   %.2f \n',...
        negElogL,kl1,kl2,kl3,klearnt1,klearnt2,klearnt3,klearnt4,kc,Prandom,km)
    
    %draw model's representations and predictions
    %--------------------------------------------
    if sum(strcmp(varargin,'vonMisesPrior'))==1
        SLdrawModelRepresentations(coh,pstd,priorModes,d,simP,PestimateGivenModel,TheModel,varargin);
    end
    
    %make sound when finished
    SLmakeSound(0.1)
    
end


%draw model predictions (centered at prior mean)
function drawMeanPreCentered(meanPred,stdPred,dataCond)

%make sure column vector
if size(meanPred,2) > size(meanPred,1)
    meanPred=meanPred';
end
if size(stdPred,2) > size(stdPred,1)
    stdPred=stdPred';
end

%factors 1,2,3
F.f1.i = dataCond(:,3);
F.f1.nm = 'd';
F.f1.L = unique(F.f1.i);
F.f1.L = sort(F.f1.L,'ascend');
F.f1.n = numel(F.f1.L);

F.f2.i = dataCond(:,2);
F.f2.nm = 'coh';
F.f2.L = unique(F.f2.i);
F.f2.L = sort(F.f2.L,'descend');
F.f2.n = numel(F.f2.L);

F.f3.i = dataCond(:,1);
F.f3.nm = 'Prior std';
F.f3.L = unique(F.f3.i);
F.f3.L = sort(F.f3.L,'descend');
F.f3.n = numel(F.f3.L);

%Graphics
F.f2.colorPre = {[0.2 0 0],...
    [0.97 0.2 0],...
    [0.8 0.4 0],...
    [0.3 0.3 0]};

%Mean data and predictions
scrsz = get(0,'ScreenSize');
fig1 = figure(1);
set(fig1,'color','w','Position',[0 1.4*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
ymax = nan(F.f2.n,1);
ymin = nan(F.f2.n,1);
h = nan(F.f2.n,1);

%graphics
ftsz = 12;
for j = 1 : F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    for i=1:F.f3.n
        
        %data for this conditions
        thisC=F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %recenter all data relative to prior
        %calculate distance to prior mean 225 deg
        xCentered = vectors2signedAngle(F.f1.i( thisC ),225);
        yCentered = vectors2signedAngle(meanPred( thisC ),225);
        
        %sort distance
        [xCenteredSorted,IA]=sort(xCentered,'ascend');
        yCenteredSorted=yCentered(IA);
        
        %sort positions and labels for plot such that 225 deg the prior
        %mean is at the center of the plot (when distance = 0).
        %Sorting is done both for x and y data such that they remain
        %aligned across conditions.
        xTickCentered = vectors2signedAngle(F.f1.L,225);
        [xTickCenteredSorted,I]=sort(xTickCentered,'ascend');
        yTickCenteredSorted = xTickCenteredSorted;
        xtickLabel=F.f1.L(I);
        ytickLabel=xtickLabel;
        
        %plot
        myPlot2 = plot(xCenteredSorted,yCenteredSorted,...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
        
        %graphics
        ymax(j,i)=max(yCenteredSorted);
        ymin(j,i)=min(yCenteredSorted);
    end
    
    %x and ylabel
    if j==1
        ylabel('Prediction of average estimate (deg)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
    
    %centered x labels
    set(gca,'ytick',yTickCenteredSorted(3:8:end),'yticklabel',ytickLabel(3:8:end))
    set(gca,'xtick',xTickCenteredSorted(3:8:end),'xticklabel',xtickLabel(3:8:end))
end

%x and ylimits
set(h,'ylim',[min(ymin(:)) max(ymax(:))])
set(h,'xlim',[min(ymin(:)) max(ymax(:))])

%clear up
SLremoveDeadSpace(0);

%selected legend
lg=legend(myPlot2);
legend(lg,'boxoff')

%----
%std
%----
fig2=figure(2);
set(fig2,'color','w','Position',[0 0 scrsz(3)/3 scrsz(4)/3])
ymax=nan(F.f2.n,1);
ymin=nan(F.f2.n,1);
h=nan(F.f2.n,1);
for j=1:F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    for i=1:F.f3.n
        
        %this conditions data
        thisC=F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %transform to distance relative to prior
        %this automatically align prior mean to the center of the plot
        %(distance=0)
        %distance to prior mean
        xCentered = vectors2signedAngle(F.f1.i( thisC ),225);
        yCentered = stdPred( thisC );
        
        %sort distance to prior mean
        [xCenteredSorted,IA]=sort(xCentered,'ascend');
        yCenteredSorted=yCentered(IA);
        
        %sort positions and labels for plot
        xTickCentered = vectors2signedAngle(F.f1.L,225);
        [xTickCenteredSorted,I]=sort(xTickCentered,'ascend');
        xtickLabel=F.f1.L(I);
        
        %plot
        myPlot2 = plot(xCenteredSorted,yCenteredSorted,...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
        
        %graphics
        xmax(j,i)=max(xCenteredSorted);
        xmin(j,i)=min(xCenteredSorted);
    end
    
    %x and ylimits
    ymin(j,i)=min(yCenteredSorted);
    
    %x and ylabel
    if j==1
        ylabel('Std of estimates (def)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
    
    %centered x labels
    set(gca,'ytick',1:4:max(stdPred),'yticklabel',1:4:max(stdPred))
    set(gca,'xtick',xTickCenteredSorted(3:8:end),'xticklabel',xtickLabel(3:8:end))
end

%x and ylimits
set(h,'ylim',[0 max(stdPred)])
set(h,'xlim',[min(xmin(:)) - 1  max(xmax(:)) + 1])

SLremoveDeadSpace(0);

%selected legend
lg=legend(myPlot2);
legend(lg,'boxoff')

%draw model predictions (not centered at prior mean)
function drawMeanPre(meanPred,stdPred,dataCond)

%make sure column vector
if size(meanPred,2) > size(meanPred,1)
    meanPred=meanPred';
end
if size(stdPred,2) > size(stdPred,1)
    stdPred=stdPred';
end

%factors 1,2,3
F.f1.i=dataCond(:,3);
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=dataCond(:,2);
F.f2.nm='coh';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=dataCond(:,1);
F.f3.nm='Prior std';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%Graphics
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};

F.f2.colorPre={[0.2 0 0],...
    [0.97 0.2 0],...
    [0.8 0.4 0],...
    [0.3 0.3 0]};

%Mean data and predictions
scrsz=get(0,'ScreenSize');
fig1=figure(1);
set(fig1,'color','w','Position',[0 1.4*scrsz(4)/3 scrsz(3)/3 scrsz(4)/3])
ymax=nan(F.f2.n,1);
ymin=nan(F.f2.n,1);
h=nan(F.f2.n,1);
count=0;

%graphics
marksz = 100;
ftsz=12;

for j=1:F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    %axis square
    for i=1:F.f3.n
        
        %this conditions data
        thisC=F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %pred
        myPlot2=plot(F.f1.i( thisC ), meanPred( thisC ),...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linestyle','-',...
            'linesmoothing','on',...
            'displayName','Bayes');
    end
    
    %x and ylimits
    ymax(j)=max(meanPred);
    ymin(j)=min(meanPred);
    xmax(j)=max(F.f1.i);
    xmin(j)=min(F.f1.i);
    
    %x and ylabel
    if j==1
        ylabel('Mean estimates (deg)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
end

%x and ylimits
set(h,'ylim',[max(ymin) max(ymax)])
set(h,'xlim',[min(xmin) max(xmax)])
SLremoveDeadSpace(0);

%selected legend
lg=legend(myPlot2);
legend(lg,'boxoff')

%-----
%std
%-----
fig2=figure(2);
set(fig2,'color','w','Position',[0 0 scrsz(3)/3 scrsz(4)/3])
ymax=nan(F.f2.n,1);
ymin=nan(F.f2.n,1);
h=nan(F.f2.n,1);

for j=1:F.f2.n
    h(j)=subplot(1,F.f2.n,j);
    %axis square
    for i=1:F.f3.n
        
        %this conditions data
        thisC=F.f2.i==F.f2.L(j) & F.f3.i==F.f3.L(i);
        hold all
        
        %pred
        myPlot2=plot(F.f1.i( thisC ), stdPred( thisC ),...
            'color',F.f2.colorPre{i},...
            'linewidth',3,...
            'linesmoothing','on',...
            'displayName','Bayes');
    end
    
    %x and ylimits
    ymax(j)=max(stdPred);
    ymin(j)=min(stdPred);
    xmax(j)=max(F.f1.i);
    xmin(j)=min(F.f1.i);
    
    %x and ylabel
    if j==1
        ylabel('Std of estimates (def)','fontsize',ftsz)
    end
    if j==round(F.f2.n/2)
        xlabel('Motion direction (deg)','fontsize',ftsz)
    end
    set(gca,'fontsize',ftsz)
end

%x and ylimits
set(h,'ylim',[max(ymin) max(ymax)])
set(h,'xlim',[min(xmin) max(xmax)])
SLremoveDeadSpace(0);

%selected legend
lg=legend(myPlot2);
legend(lg,'boxoff')


