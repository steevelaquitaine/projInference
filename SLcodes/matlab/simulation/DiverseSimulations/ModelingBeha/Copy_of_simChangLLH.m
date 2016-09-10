
%Author: Steeve Laquitaine
%date: 131024
%goal: Simulate Bayesian inference with changing LLH.
%usage: pred=simChangLLH(data,d,coh,pstd,k);
%Note: One can use the code with "drawnow" but it's 6 times faster without.

function pred=simChangLLH(data,d,coh,pstd,k)
%simulate
[pred,d,coh,pstd]=makePre(d,coh,pstd,k);
%draw
drawRawPred5(pred,data,d,coh,pstd);

%make predictions
function [pred,d,coh,pstd]=makePre(d,coh,pstd,k)
%Simulate Bayesian inference with changing LLH.
%up:  prior's mean
%kml: measurement and likelihood densities mean
%kp:  prior strength
%Pr:  probability of random estimation
%km:  motor noise strength
%P_BIweak: probability of response bias.

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
P_BIweak=k(10);

%measurement & motor noise densities
%densities used by the 3 estimation processes.
p_measStrong=vmPdfs(xe,d,kml,'norm');
p_rand=ones(numel(xe),1)/360;
mi=nan(ntrial,1);
upo1=nan(ntrial,1);
upo2=nan(ntrial,1);
upo3=nan(ntrial,1);
upo=nan(ntrial,1);
pred.sim=nan(ntrial,1);

%Calculate the percepts produced by the 3 processes.
for i=1:ntrial
    
    %1-BI with strong LIKELIHOOD
    %measurement sample -> percept
    mi(i)=randsample(xe,1,true,p_measStrong(:,i));
    upo1(i)=ra2d(de2r(mi(i),1)+atan2(sin(de2r(up,1)-de2r(mi(i),1)), (kml(i)./kp(i))'+cos(de2r(up,1)-de2r(mi(i),1))));
    
    %2-BI with weak LIKELIHOOD
    upo2(i)=225;
    
    %3-RANDOM ESTIMATES
    upo3(i)=randsample(xe,1,true,p_rand);
end

%Set mixture probabilities
%response bias trials
P_BIweaki=rand(ntrial,1)<P_BIweak;

%random estimates trials
ntrial_BIweak=sum(P_BIweaki);
ntrial_BInotweak=(ntrial-ntrial_BIweak);
Pradjusted=Pr*ntrial/ntrial_BInotweak;%works
P_BInotweaki=rand(ntrial_BInotweak,1)<Pradjusted;
Pri(P_BIweaki==0)=P_BInotweaki;

%strong LIKELIHOOD (BI) trials
P_BIstrongi=zeros(ntrial,1);
P_BIstrongi(sum([P_BIweaki Pri'],2)==0)=1;

%True values of the mixture probabilities(sum to 1)
P_BIweakReal=sum(P_BIweaki)/ntrial;
PrReal=sum(Pri)/ntrial;
P_BIstrongReal=sum(P_BIstrongi)/ntrial;

%percepts
upo(P_BIstrongi==1)=upo1(P_BIstrongi==1);
upo(P_BIweaki==1)=upo2(P_BIweaki==1);
upo(Pri==1)=upo3(Pri==1);

%estimates (add motor noise)
%calculate the estimates produced by the 3 processes.
P_motor=vmPdfs(xe,upo,km(ones(1,ntrial)),'norm');
estimate=nan(ntrial,1);
for trial=1:ntrial
    estimate(trial)=randsample(xe,1,true,P_motor(:,trial));
end

%make sure data range within 1:360.
pred.sim=estimate;
pred.sim(pred.sim==0)=360;

%calculate theoretical mean estimates
upomean=nan(ntrial,1);
for trial=1:ntrial
    mi(trial)=d(trial);
    upomean(trial)=ra2d(de2r(mi(trial),1)+atan2(sin(de2r(up,1)-de2r(mi(trial),1)), (kml(trial)./kp(trial))'+cos(de2r(up,1)-de2r(mi(trial),1))));
end

%mean is the mixture of the mean.
pred.mean=P_BIstrongReal*upomean + P_BIweakReal*up + PrReal*180;

%make sure data range within 1:360.
pred.mean(pred.mean==0)=360;

%draw
function drawRawPred5(pred,data,d,coh,pstd)
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
F.f2.color=[0.5 0 0;...
    1 0.2 0;...
    1 0.6 0;...
    0.75 0.75 0];
one=ones(size(F.f2.color,1),3);
egdec=F.f2.color+0.5*(one-F.f2.color);

%raw
figure(1);
set(gcf,'color','w')

%preallocate memory
rawData=cell(F.f1.n,F.f2.n,F.f3.n);
rawPred=cell(F.f1.n,F.f2.n,F.f3.n);
for j=1:F.f2.n
    hs(1)=subplot(2,3,j);
    hs(2)=subplot(2,3,j+3);
    for k=1:F.f1.n
        for i=1:F.f3.n
            rawData{k,i,j}=data(F.inter.pos{k,i,j},:);
            rawPred{k,i,j}=pred.sim(F.inter.pos{k,i,j},:);
            
            %raw data
            hold all
            subplot(hs(1))
            scatter(hs(1),d(F.inter.pos{k,i,j}),rawData{k,i,j}',...
                'MarkerEdgeColor',egdec(i,:),...
                'MarkerFaceColor',F.f2.color(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            alpha(0.5)
            %             drawnow
            
            %raw pred
            hold all
            subplot(hs(2))
            %colPred=F.f2.color(i,:)+0.2*([0 0 0]-F.f2.color(i,:));
            scatter(hs(2),d(F.inter.pos{k,i,j}),rawPred{k,i,j}',...
                'MarkerEdgeColor',egdec(i,:),...
                'MarkerFaceColor',F.f2.color(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
            alpha(0.5)
            %             drawnow
            
            %ideal
            subplot(hs(1))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            plot(F.f1.L,225*ones(31,1),'b.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            ylabel(hs(1),'Estimated directions (degrees)')
            
            %prior mean
            subplot(hs(2))
            plot(F.f1.L,F.f1.L,'k.','Markersize',4)
            plot(F.f1.L,225*ones(31,1),'b.','Markersize',4)
            xlim([0 360])
            ylim([0 360])
            xlabel(hs(2),'Displayed directions (degrees)')
            ylabel(hs(2),'Simulated estimates (degrees)')
        end
    end
end
title('Raw data/predictions (degrees)')
toc

%mean & std
figure(2);
set(gcf,'color','w')
zer=zeros(size(F.f2.color,1),3);
colPred=F.f2.color+0.5*(zer-F.f2.color);

%preallocate memory
statD=cell(F.f1.n,F.f2.n,F.f3.n);
statPs=cell(F.f1.n,F.f2.n,F.f3.n);
statPm=cell(F.f1.n,F.f2.n,F.f3.n);
meansD=nan(F.f1.n,F.f2.n,F.f3.n);
meansP=nan(F.f1.n,F.f2.n,F.f3.n);
stdDs=nan(F.f1.n,F.f2.n,F.f3.n);
stdPs=nan(F.f1.n,F.f2.n,F.f3.n);
for j=1:F.f2.n
    hs(1)=subplot(2,3,j);
    hs(2)=subplot(2,3,j+3);
    for i=1:F.f3.n
        for k=1:F.f1.n
            %mean & std/data & predictions
            statD{k,i,j}=circMeanStd(datac(F.inter.pos{k,i,j},:));
            meansD(k,i,j)=statD{k,i,j}.deg.mean;
            stdDs(k,i,j)=statD{k,i,j}.deg.std;
            
            %mean (theoretically)
            statPm{k,i,j}=circMeanStd(predmc(F.inter.pos{k,i,j},:));
            meansP(k,i,j)=statPm{k,i,j}.deg.mean;
            
            %std of simulation
            statPs{k,i,j}=circMeanStd(predsc(F.inter.pos{k,i,j},:));
            stdPs(k,i,j)=statPs{k,i,j}.deg.std;
        end
        
        %data mean
        hold all
        subplot(hs(1))
        hold all
        scatter(F.f1.L,meansD(:,i,j)',...
            'MarkerEdgeColor',egdec(i,:),...
            'MarkerFaceColor',F.f2.color(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        drawnow
        
        %pred mean
        plot(F.f1.L,meansP(:,i,j)','-',...
            'linewidth',2,...
            'color',colPred(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        box off
        %         drawnow
        
        %data std
        hold all
        subplot(hs(2))
        hold all
        scatter(F.f1.L,stdDs(:,i,j)',...
            'MarkerEdgeColor',egdec(i,:),...
            'MarkerFaceColor',F.f2.color(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        %         drawnow
        
        %pred std
        plot(F.f1.L,stdPs(:,i,j)','-',...
            'linewidth',2,...
            'color',colPred(i,:),...
            'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
        box off
        %         drawnow
        
        %ideal
        subplot(hs(1))
        plot(F.f1.L,F.f1.L,'k.','Markersize',4)
        plot(F.f1.L,225*ones(31,1),'b.','Markersize',4)
        xlim([0 360])
        ylim([0 360])
        ylabel(hs(1),'Estimated directions (degrees)')
        %         drawPublishAxis
        
        %std
        subplot(hs(2))
        plot(225*ones(31,1),F.f1.L,'b.','Markersize',4)
        xlim([0 360])
        ylim([0 200])
        xlabel(hs(2),'Displayed directions (degrees)')
        ylabel(hs(2),'Std of estimated directions (degrees)')
        %         drawPublishAxis
    end
end
title('Mean/std of data/predictions (degrees)')

%support
%von mises
function mPdfs=vmPdfs(x,u,k,type)
%support
%Von Mises
%author: Steeve Laquitaine
%usage:
%mPdfs=vmPdfs(x,u,k,'norm')

%inputs:
%'norm' scales to probabilities, otherwise [];

%Description: von Mises
%This code produces von Mises densities vm(u,k) based on the equation
%vm=exp(k.*cos(x-u))./(2*pi.*besseli(0,k));
%The equation is adjusted because of the following numerical issues:
%when k>700, vm is NaN because besseli(0,k) and exp(k.*cos(x-u)) reach numerical limits.
%exp(k.*cos(x-u)-k) scales vm without changing its shape.
%besseli(0,k,1)) does same.
%The adjusted equation and the exact equation yield exact same results
%except that the adjusted equation works for large k.
w=1;

%radians
x=de2r(x,1); x=x';
u=de2r(u,1); u=u';
x2=x(:,ones(numel(u),1));
u2=u(ones(numel(x),1),:);
k2=k(ones(numel(x),1),:);

%von mises
mPdfs=exp(k2.*cos(w*(x2-u2))-k2)./(2*pi.*besseli(0,k2,1));

%scale to pdfs
if strcmp(type,'norm')==1
    Z_=sum(mPdfs);
    Z=Z_(ones(numel(x),1),:);
    mPdfs=mPdfs./Z;
else
end