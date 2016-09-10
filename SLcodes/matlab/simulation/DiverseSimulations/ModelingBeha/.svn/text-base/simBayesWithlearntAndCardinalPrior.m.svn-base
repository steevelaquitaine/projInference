
%Author: Steeve Laquitaine
%date  : 140205
%goal  : Simulate Bayesian inference with changing LLH.
%usage :[pred,Kinput,Kreal,d,coh,pstd]=simBayesWithlearntAndCardinalPrior(data,d,coh,pstd,k,'TheoreticalMean',10)
%[pred,Kinput,Kreal,d,coh,pstd]=simBayesWithlearntAndCardinalPrior(data,d,coh,pstd,k,'TrueMean',10)
%Note: One can use the code with "drawnow" but it's 6 times faster without.
%description:

function [pred,Kinput,Kreal,d,coh,pstd]=simBayesWithlearntAndCardinalPrior(data,d,coh,pstd,k,meanType,numsim)

%get more samples by repeating simulations 'numsim' times.
data=repmat(data,numsim,1);
d=repmat(d,numsim,1);
coh=repmat(coh,numsim,1);
pstd=repmat(pstd,numsim,1);

%simulate
[pred,d,coh,pstd,Kinput,Kreal]=makePre(d,coh,pstd,k);
fprintf('calculating predictions.')

%draw
drawRawPre(pred,data,d,coh,pstd,meanType,numsim);
fprintf('drawing ...')

%make predictions
function [pred,d,coh,pstd,Kinput,Kreal]=makePre(d,coh,pstd,k)
%Simulate Bayesian inference with cardinal and learnt priors
%kml    :measurement and likelihood densities strengths
%kcard  :cardinal prior strength
%Prandom:probability of random estimation
%kmotor :motor noise strength
%PlearntPriorSwitch: probability of response bias.

%useful variables
ntrial=numel(d);
xe=1:1:360;

%strengths of measurement distributions, learnt, cardinal prior, motor 
%noise; probability of random estimation, and of switch to learnt prior
kml(coh==0.24)    =k(1);
kml(coh==0.12)    =k(2);
kml(coh==0.06)    =k(3);
klearnt(pstd==80) =k(4);
klearnt(pstd==40) =k(5);
klearnt(pstd==20) =k(6);
klearnt(pstd==10) =k(7);
kcard(1:ntrial)   =k(8);
Prandom           =k(9);
kmotor            =k(10);
PlearntPriorSwitch=k(11);

%measurement & motor noise densities
%densities used by the 3 estimation processes.
Pmeas=vmPdfs(xe,d,kml,'norm');
p_rand=ones(numel(xe),1)/360;
upo1=nan(ntrial,1);
upo2=nan(ntrial,1);
upo3=nan(ntrial,1);
upo=nan(ntrial,1);
pred.sim=nan(ntrial,1);
kcardi=kcard(1);
up=225;

%cardinal prior
vmc=vmPdfs(xe,[90 180 270 360],kcardi,'norm');
PRIORcardinal=0.25.*sum(vmc,2);
PRIORcardinal=PRIORcardinal(:,ones(ntrial,1));

%learnt prior
%PRIORlearnt=vmPdfs(xe,up(ones(ntrial),1),klearnt,'norm');
PRIORlearnt=vmPdfs(xe,up(ones(1,ntrial)),klearnt,'norm');



%Calculate the percepts produced by the two processes: Bayesian inference
%and random estimation.
%Parfor is very useful here. It speeds the code by a factor of 7.
parfor i=1:ntrial
    
    %1-Bayesian inference
    %measurement sample -> percept
    mi=randsample(xe,1,true,Pmeas(:,i));
    llh=vmPdfs(xe,mi,kml(i),'norm');
    
    %posterior: probability of the common causes of llh and priors
    posterior=llh.*PRIORcardinal(:,i).*PRIORlearnt(:,i)./sum(llh.*PRIORcardinal(:,i).*PRIORlearnt(:,i));
    
    %calculate percept upo, the most likely common cause.
    %How to choose 'upo' when posterior has many modes?
    %The posterior has two modes at max and it happens when motion
    %direction is at oblique directions. We assume that at the obliques
    %subjects estimates the two modes of the posterior and randomly choose
    %between them. When llh is flat the modes are the nearest cardinals
    %and when llh is non-flat the modes are in between the oblique and the
    %nearest cardinals.
    %I could not find two modes at the obliques because of resolution limits.
    %(number of floating points). So I reasoned that the one mode that is 
    %not found can be identified knowing that it has to be at equal 
    %distance to the oblique. But adding learnt prior complicate more the
    %predictions. e.g., now when learnt prior equal llh strength and
    %llh mean is the opposite of the prior mean (e.g., 45 degrees), llh and
    %learnt prior cancel each other and posterior peaks at cardinal modes.
    %So this doesn't hold anymore and my not hold for other situation i
    %cannot think of..so I rather try with resolution low enough to include
    %the lowest reasonable resolution possible (when llh and priors are 
    %very weak).
    posterior=round(posterior.*10^6)./10^6;
    MAPs=find(posterior==max(posterior));
    equiproba=1/numel(MAPs)*ones(numel(MAPs),1);
    
    %We add zero because 'randsample' does not work with single-value MAPs
    upo1(i)=randsample([MAPs;0],1,true,[equiproba;0]); 
end

%2-switch to prior
upo2=up(ones(ntrial,1));

%3-random estimation
upo3=randsample(xe,ntrial,true,p_rand)';

%Set mixture probabilities
%Switch to prior trials. When set to 0 there is no switch.
PlearntPriorSwitchi=rand(ntrial,1)<PlearntPriorSwitch;

%random estimates trials
Prandomi=zeros(ntrial,1);
ntrial_BIweak=sum(PlearntPriorSwitchi);
ntrial_BInotweak=(ntrial-ntrial_BIweak);
Prandomadjusted=Prandom*ntrial/ntrial_BInotweak;%works
P_BInotweaki=rand(ntrial_BInotweak,1)<Prandomadjusted;
Prandomi(PlearntPriorSwitchi==0)=P_BInotweaki;

%Bayesian trials
P_BIstrongi=zeros(ntrial,1);
P_BIstrongi(sum([PlearntPriorSwitchi Prandomi],2)==0)=1;

%true values of the mixture probabilities(sum to 1)
mixtP.PlearntPriorSwitchReal=sum(PlearntPriorSwitchi)/ntrial;
mixtP.PrandomReal=sum(Prandomi)/ntrial;
mixtP.P_BIstrongReal=sum(P_BIstrongi)/ntrial;

%percepts
upo(P_BIstrongi==1)=upo1(P_BIstrongi==1);
upo(PlearntPriorSwitchi==1)=upo2(PlearntPriorSwitchi==1);
upo(Prandomi==1)=upo3(Prandomi==1);

%estimates (add motor noise)
%calculate the estimates produced by the 3 processes.
P_motor=vmPdfs(xe,upo,kmotor(ones(1,ntrial)),'norm');
estimate=nan(ntrial,1);
parfor trial=1:ntrial
    estimate(trial)=randsample(xe,1,true,P_motor(:,trial));
end

%make sure data range within 1:360.
pred.sim=estimate;
pred.sim(pred.sim==0)=360;

%calculate theoretical mean estimates over trial (Bayesian mean estimates).
%Motor noise is ignored because it averages out over trials and mean
%estimates peak at the value of the percept. Random estimation is also
%ignored because the mean of a uniform circular density is null vector.
%The mean estimate is the vector average over of the posterior's mode(s). 
%We assumed that at the obliques subjects randomly sample the posterior's 
%modes.
%We checked that the code run fine when we round posterior at a precision 
%of 10^-6 by set flat learnt priors and by looking at if the mean estimate 
%at the obliques was the obliques. 
predmean=nan(ntrial,1);

%llh
llh=vmPdfs(xe,d,kml,'norm');

%posterior
posterior=llh.*PRIORcardinal.*PRIORlearnt;
Z=sum(posterior,1);
Z=Z(ones(numel(xe),1),:);
posterior=posterior./Z;
posterior=round(posterior.*10.^6)./10.^6;
for trial=1:ntrial
    MAPs=find(posterior(:,trial)==max(posterior(:,trial)));
    predmean(trial)=vectMean(MAPs,ones(1,numel(MAPs)));
end
pred.mean=predmean;

%make sure data range within 1:360.
pred.mean(pred.mean==0)=360;

%Re-adjust the fit parameters with the true mixture probabilities and not
%the input.
%backup input
Kinput=k;

%adjust Probability of random estimations "Prandom" and probability of
%switching to the prior "PlearntPriorSwitch".
Kreal=k;
Kreal(9)=mixtP.PrandomReal;
Kreal(11)=mixtP.PlearntPriorSwitchReal;

%draw
function drawRawPre(pred,data,d,coh,pstd,meanType,numsim)
tic
%convert from polar to cartesian coordinates to calculate the statistics of
%subject's data and of the statistics of the simulation.
datac=polar2cartesian(data,1);
predsc=polar2cartesian(pred.sim,1);
predmc=polar2cartesian(pred.mean,1);

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
egdec=F.f2.color0+0.8*(one-F.f2.color0);

%predictions
zer=zeros(size(F.f2.color0,1),3);
colPrandomed=F.f2.color0+0*(zer-F.f2.color0);
F.f2.color=F.f2.color0+0*(one-F.f2.color0);
markersizeS=50;
%line
linWdth=4;
%font
ft=14;

%RAW
figure(1);
set(gcf,'color','w')

%plot only one instance of the simulations with the raw data (for clarity).
pred.sim1=nan(numel(d),1);
pred.sim1(1:numel(d)/numsim)=pred.sim(1:numel(d)/numsim);


Finterpos=F.inter.pos;
tic
for j=1:F.f2.n
    hs(1)=subplot(2,3,j);
    hs(2)=subplot(2,3,j+3);
    for k=1:F.f1.n
        for i=1:F.f3.n
            pos=Finterpos{k,i,j};
            rawData=data(pos,:);
            rawPrandomed=pred.sim1(pos,:);
            
            %raw data
            hold all
            subplot(hs(1))
            scatter(hs(1),d(pos),rawData',...
                'MarkerEdgeColor',egdec(i,:),...
                'MarkerFaceColor',F.f2.color0(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            
            
            %raw pred
            hold all
            subplot(hs(2))
            scatter(hs(2),d(pos),rawPrandomed',...
                'MarkerEdgeColor',egdec(i,:),...
                'MarkerFaceColor',F.f2.color0(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            
            %ideal
            subplot(hs(1))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            plot(F.f1.L,225*ones(numel(F.f1.L),1),'b.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            ylabel(hs(1),'Estimated directions (degrees)')
            
            subplot(hs(2))
            %ideal
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            
            %prior mean
            plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4)
            
            %cardinals
            plot(F.f1.L,90*ones(F.f1.n,1),'k.','Markersize',4)
            plot(F.f1.L,180*ones(F.f1.n,1),'k.','Markersize',4)
            plot(F.f1.L,270*ones(F.f1.n,1),'k.','Markersize',4)
            plot(F.f1.L,360*ones(F.f1.n,1),'k.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            xlabel(hs(2),'Displayed directions (degrees)')
            ylabel(hs(2),'Simulated estimates (degrees)')
        end
    end
end
toc


%mean & std
figure(2);
set(gcf,'color','w')

%preallocate memory
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
            %mean & std of data
            statD{k,i,j}=circMeanStd(datac(F.inter.pos{k,i,j},:));
            meansD(k,i,j)=statD{k,i,j}.deg.mean;
            stdDs(k,i,j)=statD{k,i,j}.deg.std;
            
            %mean & std of simulation
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
            'MarkerEdgeColor',egdec(i,:),...
            'MarkerFaceColor',F.f2.color(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        %Prandomedictions' mean
        if strcmp(meanType,'TheoreticalMean')
            %theoretical
            plot(F.f1.L,meansP(:,i,j)','-',...
                'linewidth',linWdth,...
                'color',colPrandomed(i,:),...
                'lineSmoothing','on',...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        elseif strcmp(meanType,'TrueMean')
            %true
            plot(F.f1.L,meansPs(:,i,j)','-',...
                'linewidth',linWdth,...
                'color',colPrandomed(i,:),...
                'lineSmoothing','on',...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        end
        
        %data std
        hold all
        subplot(hs(2))
        hold all
        scatter(F.f1.L,stdDs(:,i,j)',...
            markersizeS,...
            'MarkerEdgeColor',egdec(i,:),...
            'MarkerFaceColor',F.f2.color(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        %predicted std
        plot(F.f1.L,stdPs(:,i,j)','-',...
            'linewidth',linWdth,...
            'color',colPrandomed(i,:),...
            'lineSmoothing','on',...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        
        subplot(hs(1))
        %ideal predictions
        plot(F.f1.L,F.f1.L,'k.','Markersize',4,'lineSmoothing','on')
        
        %prior mean
        plot(F.f1.L,225*ones(F.f1.n,1),'b.','Markersize',4,'lineSmoothing','on')
        
        %cardinals
        plot(F.f1.L,90*ones(F.f1.n,1),'k.','Markersize',4,'lineSmoothing','on')
        plot(F.f1.L,180*ones(F.f1.n,1),'k.','Markersize',4,'lineSmoothing','on')
        plot(F.f1.L,270*ones(F.f1.n,1),'k.','Markersize',4,'lineSmoothing','on')
        plot(F.f1.L,360*ones(F.f1.n,1),'k.','Markersize',4,'lineSmoothing','on')
        xlim([30 365])
        ylim([30 365])
        set(gca,'xtick',[30 225 365],'xticklabel',{'30','225','365'},'fontsize',ft)
        set(gca,'ytick',[30 225 365],'yticklabel',{'30','225','365'},'fontsize',ft)
        ylabel(subplot(2,3,1),'Mean estimated directions (degrees)')
        
        %std
        subplot(hs(2))
        plot(225*ones(F.f1.n,1),F.f1.L,'b.','Markersize',5,'lineSmoothing','on')
        xlim([30 365])
        ylim([0 150])
        set(gca,'xtick',[30 225 365],'xticklabel',{'30','225','365'},'fontsize',ft)
        xlabel(subplot(2,3,5),'Displayed directions (degrees)')
        ylabel(subplot(2,3,4),'Std of estimated directions (degrees)')
    end
end
