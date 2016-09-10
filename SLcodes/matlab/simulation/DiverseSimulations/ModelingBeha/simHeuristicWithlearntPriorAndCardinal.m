
%Author: Steeve Laquitaine
  %date: 140124
  %goal: Simulate "Heuristic" model where likelihood mean and prior mean are
  %biased toward cardinal prior.
 %usage: [pred,Kinput,Kreal,mixProba]=simHeuristic(data,d,coh,pstd,k,'TrueMean');
  %Note: One can use the code with "drawnow" but it's 6 times faster without.
  
function [pred,Kinput,Kreal,mixProba]=simHeuristicWithlearntPriorAndCardinal(data,d,coh,pstd,k,meanType,numsim)
%repeat simulations 'numsim' times
data=repmat(data,numsim,1);
d=repmat(d,numsim,1);
coh=repmat(coh,numsim,1);
pstd=repmat(pstd,numsim,1);

%simulate
[pred,d,coh,pstd,Kinput,Kreal,mixProba]=makePre(d,coh,pstd,k);
fprintf('calculating predictions.')

%draw
drawRawPred5(pred,data,d,coh,pstd,meanType,numsim);
fprintf('drawing ...')

%make predictions
function [pred,d,coh,pstd,Kinput,Kreal,mixProba]=makePre(d,coh,pstd,k)
%Simulate Bayesian inference with heuristic model.
%up:  prior's mean
%kml: measurement and likelihood densities mean
%kp:  prior strength
%Pr:  probability of random estimation
%km:  motor noise strength
%Ppriormn: probability of response bias.

%useful variables
ntrial=numel(d);
xe=1:1:360;

%strength of the measurement distributions, priors & motor noise
up=225;
kml(coh==0.24)=k(1);
kml(coh==0.12)=k(2);
kml(coh==0.06)=k(3);
kp(pstd==80)=k(4);
kp(pstd==40)=k(5);
kp(pstd==20)=k(6);
kp(pstd==10)=k(7);
Pr=k(8);
km=k(9);
Beta=k(10);
kcard=k(11);

%measurement & motor noise densities
%densities used by the 3 estimation processes.
p_meas=vmPdfs(xe,d,kml,'norm');
p_rand=ones(numel(xe),1)/360;
mi=nan(ntrial,1);
upo1=nan(ntrial,1);
upo2=nan(ntrial,1);
upo3=nan(ntrial,1);
upo=nan(ntrial,1);
pred.sim=nan(ntrial,1);


%Calculate the percepts produced by the 3 processes.
%Cardinal prior
vm=vmPdfs(xe,[90 180 270 360],kcard,'norm');
priorCardinal=0.25.*(vm(:,1)+vm(:,2)+vm(:,3)+vm(:,4));
priorCardinal=priorCardinal(:,ones(ntrial,1));

%learnt
priorLearnt=vmPdfs(xe,up(ones(ntrial,1)),kp,'norm');
% priorLearnt=priorLearnt(:,ones(ntrial,1));

for i=1:ntrial
    
    %1. perceive likelihood mean ('seeing') biased by cardinals (BI)
    %The posterior has two modes at max and it happens when motion
    %direction is at oblique directions. We assume that at the obliques
    %subjects estimates the two modes of the posterior and random choose
    %between them.
    mi(i)=randsample(xe,1,true,p_meas(:,i));    
    llh(:,i)=vmPdfs(xe,mi(i),kml(i),'norm');
    posteriorforllh(:,i)=llh(:,i).*priorCardinal(:,i)/sum(llh(:,i).*priorCardinal(:,i));
    posteriorforllh(:,i)=fix(posteriorforllh(:,i)*10^10)/10^10;
    
    %find MAPs estimates (values) of each ul/mi (rows)
    numupo1=sum(posteriorforllh(:,i)==max(posteriorforllh(:,i)));
    thisupo1(1:numupo1)=xe(posteriorforllh(:,i)==max(posteriorforllh(:,i)));
    randsampleupo1=randperm(numel(thisupo1));
    upo1(i)=thisupo1(randsampleupo1(1));
    
    %2. perceive learnt prior mean biased by cardinals (BI)
    posteriorforlearntPrior(:,i)=priorCardinal(:,i).*priorLearnt(:,i)/sum(priorLearnt(:,i).*priorCardinal(:,i));
    posteriorforlearntPrior(:,i)=fix(posteriorforlearntPrior(:,i)*10^10)/10^10;
    
    %find MAPs estimates (values) of each ul/mi (rows)
    numupo2=sum(posteriorforlearntPrior(:,i)==max(posteriorforlearntPrior(:,i)));
    thisupo2(1:numupo2)=xe(posteriorforlearntPrior(:,i)==max(posteriorforlearntPrior(:,i)));
    randsampleupo2=randperm(numel(thisupo2));
    upo2(i)=thisupo2(randsampleupo2(1));

    %3. make random guess
    upo3(i)=randsample(xe,1,true,p_rand);
end


%Compute fractions of trials controlled by the three processes. The 
%probability that percepts are attracted toward prior increases as the 
%ratio of prior and likelihood strengths increase and the probability of
%random estimation is fixed and low. We define them as "weights" first.
%when Beta is high (e.g., >=1 denominator is inf for large kml and kp..) 
%and there are other cases like that where the equation suffers numerical
%problems. So I scale kp and kml with log which basically is getting rid of
%exp.
%weightPriormn=exp(log(Beta.*kp))./(exp(log(Beta.*kp))+exp(log(Beta.*kml)));
%--> weightPriormn=Beta.*kp./(Beta.*(kp+kml));
%But if we do that we lose Beta too that we need to control sensitivity to
%prior and llh strength. So we use:
%kp-independent
weightPriormn=exp(Beta.*log(kp))./(exp(Beta.*log(kp))+exp(Beta.*log(kml)));
weightLlhmn=1-weightPriormn;

%scale "weights" to probabilities
sumP=weightPriormn+weightLlhmn+Pr;
Ppriormn=1-(Pr+weightLlhmn)./sumP;
Pllhmn=1-(Pr+weightPriormn)./sumP;

%assign each process to a trial.
%trials biased to prior mean
Ppriormni=rand(ntrial,1)<Ppriormn';

%trials with random estimates
Pri=zeros(ntrial,1);
ntrial_prior=sum(Ppriormni);
ntrial_llh=(ntrial-ntrial_prior);
Pradjusted=Pr*ntrial/ntrial_llh;
Pritmp=rand(ntrial_llh,1)<Pradjusted;
Pri(Ppriormni==0)=Pritmp;

%trials biased to likelihood mean
Pllhmni=zeros(ntrial,1);
Pllhmni(sum([Ppriormni Pri],2)==0)=1;

%true mixture probabilities(sum to 1)
%get the probabilities that each process is in control for each of the 3x4
%conditions (coh x prior)
cohVal=unique(coh);
mixProba.PpriormnReal=nan(4,numel(cohVal));
mixProba.PllhmnReal=nan(4,numel(cohVal));
mixProba.colname={'coh0.06','coh 0.12','coh 0.24'};
mixProba.rwname={'pstd10','pstd20','pstd40','pstd80'};
for i=1:numel(cohVal)
    cohi=coh==cohVal(i);
    
    %bias to prior
    mixProba.PpriormnReal(1,i)=sum(Ppriormni(cohi&pstd==10))/sum(cohi&pstd==10);
    mixProba.PpriormnReal(2,i)=sum(Ppriormni(cohi&pstd==20))/sum(cohi&pstd==20);
    mixProba.PpriormnReal(3,i)=sum(Ppriormni(cohi&pstd==40))/sum(cohi&pstd==40);
    mixProba.PpriormnReal(4,i)=sum(Ppriormni(cohi&pstd==80))/sum(cohi&pstd==80);
    
    %bias to llh
    mixProba.PllhmnReal(1,i)=sum(Pllhmni(cohi&pstd==10))/sum(cohi&pstd==10);
    mixProba.PllhmnReal(2,i)=sum(Pllhmni(cohi&pstd==20))/sum(cohi&pstd==20);
    mixProba.PllhmnReal(3,i)=sum(Pllhmni(cohi&pstd==40))/sum(cohi&pstd==40);
    mixProba.PllhmnReal(4,i)=sum(Pllhmni(cohi&pstd==80))/sum(cohi&pstd==80);
    
    %random
    mixProba.PrReal(1,i)=sum(Pri(cohi&pstd==10))/sum(cohi&pstd==10);
    mixProba.PrReal(2,i)=sum(Pri(cohi&pstd==20))/sum(cohi&pstd==20);
    mixProba.PrReal(3,i)=sum(Pri(cohi&pstd==40))/sum(cohi&pstd==40);
    mixProba.PrReal(4,i)=sum(Pri(cohi&pstd==80))/sum(cohi&pstd==80);
end
Prreal=sum(Pri)/ntrial;

%time series of percepts
upo(Pllhmni==1)=upo1(Pllhmni==1);
upo(Ppriormni==1)=upo2(Ppriormni==1);
upo(Pri==1)=upo3(Pri==1);

%estimates (add motor noise)
%calculate the estimates produced by the 3 processes.
P_motor=vmPdfs(xe,upo,km(ones(1,ntrial)),'norm');
estimate=nan(ntrial,1);
parfor trial=1:ntrial
    estimate(trial)=randsample(xe,1,true,P_motor(:,trial));
end

%make sure data range within 1:360





!!!!!!!!!


pred.sim=estimate;
pred.sim(pred.sim==0)=360;

%calculate theoretical mean estimates over trial.
%mean estimate is the mixture of the mean estimates of the 3 processes.
%Motor noise is ignored because it averages out over trials and mean
%estimates peak at the value of the percept.
upomnllh=nan(ntrial,1);
predmn=nan(ntrial,1);
for trial=1:ntrial
    %percept=measurement sample
    mi(itrial)=d(trial);    
    llh(:,trial)=vmPdfs(xe,mi(trial),kml(trial),'norm');
    posteriorforllh(:,trial)=llh(:,trial).*priorCardinal(:,itrial)/sum(llh(:,trial).*priorCardinal(:,trial));
    posteriorforllh(:,trial)=fix(posteriorforllh(:,trial)*10^10)/10^10;
        
    %find MAPs estimates (values) of each ul/mi (rows)
    numupo1=sum(posteriorforllh(:,trial)==max(posteriorforllh(:,trial)));
    thisupo1(1:numupo1)=xe(posteriorforllh(:,trial)==max(posteriorforllh(:,trial)));
    upomnllh(trial)=mean(thisupo1);
    
    %model's mn estimate
    %Here we need to ignore the trials when choice are random because they
    %don't affect the mean and it the mean of a uniform von mises 
    %distribution is not clear/unknown.
    %When we switch to the prior biased by cardinal prior, the average
    %estimate is still the prior mean.
    predmn(trial)=vectMean([upomnllh(trial);up],[Pllhmn(trial),Ppriormn(trial)]);
end
pred.mn=predmn;

%make sure data range within 1:360.
pred.mn(pred.mn==0)=360;

%Re-adjust the fit parameters with the true mixture probabilities and not
%the ones input.
%backup input
Kinput=k;

%adjust Probability of random estimations "Pr" and probability of weak
%likelihood "Ppriormn".
Kreal=k;
Kreal(8)=Prreal;
Kreal(10)=Beta;

%draw
function drawRawPred5(pred,data,d,coh,pstd,meanType,numsim)
%convert from polar to cartesian coordinates to calculate the statistics of
%subject's data and of the statistics of the simulation.
datac=polar2cartesian(data,1);
predsc=polar2cartesian(pred.sim,1);
predmc=polar2cartesian(pred.mn,1);

%factors 1,2,3
F.f1.i=d;
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=coh;
F.f2.nm='coh';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=pstd;
F.f3.nm='pstd';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%positions main
for i=1:F.f1.n
    F.f1.pos(i)={find(F.f1.i==F.f1.L(i))};
end
for i=1:F.f2.n
    F.f2.pos(i)={find(F.f2.i==F.f2.L(i))};
end
for i=1:F.f3.n
    F.f3.pos(i)={find(F.f3.i==F.f3.L(i))};
end

%positions inter
for k=1:F.f1.n
    for j=1:F.f2.n
        for i=1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end

%Graphics
F.f2.color0=[0.5 0 0;...
    1 0.2 0;...
    1 0.6 0;...
    0.75 0.75 0];
one=ones(size(F.f2.color0,1),3);

%marker's edges, color & size
egdec=F.f2.color0+1*(one-F.f2.color0);
egdecrw=F.f2.color0+0.3*(one-F.f2.color0);
zer=zeros(size(F.f2.color0,1),3);
colPred=F.f2.color0+0.3*(zer-F.f2.color0);
F.f2.color=F.f2.color0+0.5*(one-F.f2.color0);
markersizeS=40;

%line
linWdth=3;

%font
ft=14;

%RAW data & simulations
figure;
set(gcf,'color','w')
rawData=cell(F.f1.n,F.f2.n,F.f3.n);
rawPred=cell(F.f1.n,F.f2.n,F.f3.n);

%plot only one instance of the simulations with the raw data (for clarity).
pred.sim1=nan(numel(d),1);
pred.sim1(1:numel(d)/numsim)=pred.sim(1:numel(d)/numsim);
for j=1:F.f2.n
    hs(1)=subplot(2,3,j);
    hs(2)=subplot(2,3,j+3);
    for k=1:F.f1.n
        for i=1:F.f3.n
            rawData{k,i,j}=data(F.inter.pos{k,i,j},:);
            rawPred{k,i,j}=pred.sim1(F.inter.pos{k,i,j},:);
            
            %data
            hold all
            subplot(hs(1))
            scatter(hs(1),d(F.inter.pos{k,i,j}),rawData{k,i,j}',...
                'MarkerEdgeColor',egdecrw(i,:),...
                'MarkerFaceColor',F.f2.color0(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            
            %simulations
            hold all
            subplot(hs(2))
            %colPred=F.f2.color(i,:)+0.2*([0 0 0]-F.f2.color(i,:));
            scatter(hs(2),d(F.inter.pos{k,i,j}),rawPred{k,i,j}',...
                'MarkerEdgeColor',egdecrw(i,:),...
                'MarkerFaceColor',F.f2.color0(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            
            %ideal
            subplot(hs(1))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            ylabel(hs(1),{'Estimated directions', '(degrees)'},'fontsize',ft)
            
            %prior mean
            subplot(hs(2))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            xlabel(hs(2),{'Displayed directions','(degrees)'},'fontsize',ft)
            ylabel(hs(2),{'Simulated estimates','(degrees)'},'fontsize',ft)
            %drawPublishAxis
        end
    end
end
toc


%MEAN & STD of data & simulations
figure;
set(gcf,'color','w')
statD=cell(F.f1.n,F.f2.n,F.f3.n);
statPs=cell(F.f1.n,F.f2.n,F.f3.n);
statPm=cell(F.f1.n,F.f2.n,F.f3.n);
meansD=nan(F.f1.n,F.f2.n,F.f3.n);
meansP=nan(F.f1.n,F.f2.n,F.f3.n);
meansPs=nan(F.f1.n,F.f2.n,F.f3.n);
stdDs=nan(F.f1.n,F.f2.n,F.f3.n);
stdPs=nan(F.f1.n,F.f2.n,F.f3.n);
for j=1:F.f2.n
    hs(1)=subplot(2,3,j);
    hs(2)=subplot(2,3,j+3);
    for i=1:F.f3.n
        for k=1:F.f1.n
            %data
            statD{k,i,j}=circMeanStd(datac(F.inter.pos{k,i,j},:));
            meansD(k,i,j)=statD{k,i,j}.deg.mean;
            stdDs(k,i,j)=statD{k,i,j}.deg.std;
            
            %simulation
            statPs{k,i,j}=circMeanStd(predsc(F.inter.pos{k,i,j},:));
            meansPs(k,i,j)=statPs{k,i,j}.deg.mean;
            stdPs(k,i,j)=statPs{k,i,j}.deg.std;
            
            %mean (theoretically)
            statPm{k,i,j}=circMeanStd(predmc(F.inter.pos{k,i,j},:));
            meansP(k,i,j)=statPm{k,i,j}.deg.mean;
        end
        
        %data mean
        hold all
        subplot(hs(1))
        hold all
        scatter(F.f1.L,meansD(:,i,j)',...
            markersizeS,...
            'MarkerEdgeColor',egdecrw(i,:),...
            'MarkerFaceColor',F.f2.color0(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        %Predictions' mean
        if strcmp(meanType,'TheoreticalMean')
            %theoretical
            plot(F.f1.L,meansP(:,i,j)','-',...
                'linewidth',linWdth,...
                'color',colPred(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        elseif strcmp(meanType,'TrueMean')
            %true
            plot(F.f1.L,meansPs(:,i,j)','-',...
                'linewidth',linWdth,...
                'color',colPred(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        end
        
        %data std
        hold all
        subplot(hs(2))
        hold all
        scatter(F.f1.L,stdDs(:,i,j)',...
            markersizeS,...
            'MarkerEdgeColor',egdecrw(i,:),...
            'MarkerFaceColor',F.f2.color0(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        %predicted std
        plot(F.f1.L,stdPs(:,i,j)','-',...
            'linewidth',linWdth,...
            'color',colPred(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        %ideal predictions
        subplot(hs(1))
        plot(F.f1.L,F.f1.L,'k.','Markersize',4,'lineSmoothing','on')
        plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4,'lineSmoothing','on')
        xlim([30 365])
        ylim([30 365])
        set(gca,'xtick',[30 225 365],'xticklabel',{'30','225','365'},'fontsize',ft)
        set(gca,'ytick',[30 225 365],'yticklabel',{'30','225','365'},'fontsize',ft)
        ylabel(subplot(2,3,1),{'Mean estimated directions','(degrees)'})
        
        %std
        subplot(hs(2))
        plot(225*ones(F.f1.n,1),F.f1.L,'b.','Markersize',5,'lineSmoothing','on')
        xlim([30 365])
        ylim([0 150])
        set(gca,'xtick',[30 225 365],'xticklabel',{'30','225','365'},'fontsize',ft)
        xlabel(subplot(2,3,5),{'Displayed directions','(degrees)'})
        ylabel(subplot(2,3,4),{'Std of estimated directions','(degrees)'})
    end
end
