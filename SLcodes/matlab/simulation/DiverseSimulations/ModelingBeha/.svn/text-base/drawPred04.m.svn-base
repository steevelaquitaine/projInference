
%Author: steeve Laquitaine
%date: 131024

%usage: pred=drawRawDataFitwithKm(data,disp,coh,pstd,k,numel(data));

%Note: One can use the code with "drawnow" but it's 6 times faster without.


function pred=drawPred04(data,d,coh,pstd,k)
pred=makePre(d,coh,pstd,k);
drawRawPred5(pred,data,d,coh,pstd);

%Make predictions
function pred=makePre(d,coh,pstd,k)

%widths of the measurement distributions, priors & motor noise
up=225;
kml(coh==0.24)=k(1);
kml(coh==0.12)=k(2);
kml(coh==0.06)=k(3);
kp(pstd==80)=k(4);
kp(pstd==40)=k(5);
kp(pstd==20)=k(6);
kp(pstd==10)=k(7);
Pm=k(8);
km=k(9);
Pp=k(10);

%measurement & motor noise densities
xe=1:1:360;

%Calculate the probability density of producing data at each trial
%(same as the llh of observing data given the model) and sample from this
%density to simulate data. Here the mean of llh is d and not mi because it
%d is the most common measurement (mean of the measurement distribution). 
%that should be sampled in a simulation.
%llh
llh=vmPdfs(xe,d,kml,'norm');

%or response bias
prior=vmPdfs(xe,up(ones(numel(d),1)),kp,'norm');
lorRB=(1-Pp).*llh+Pp.*prior;

%and motor noise
Pmot=vmPdfs(xe,0,km,'norm');
Pmot=Pmot(:,ones(size(lorRB,2),1));
ML.e=circConv(lorRB,Pmot);

%or random choice
Pe=(1-Pm).*ML.e+Pm;

%scale to probability
Z_=sum(Pe);
Z=Z_(ones(size(Pe,1),1),:);
Pe=abs(Pe./Z);

%preallocate
mi=nan(numel(d),1);
upo=nan(numel(d),1);
pred.sim=nan(numel(d),1);
pred.mean=nan(numel(d),1);
for i=1:numel(d)
    
    %trial-estimates
    pred.sim(i)=randsample(xe,1,true,Pe(:,i));
  
    %trial-average (theoretically)
    %perceived
    mi(i)=d(i);
    upo(i)=ra2d(de2r(mi(i),1)+atan2(sin(de2r(up,1)-de2r(mi(i),1)), (kml(i)./kp(i))'+cos(de2r(up,1)-de2r(mi(i),1))));
    %with response bias
    pred.mean(i)=(1-Pp).*upo(i)+Pp.*up;
end

%make sure data range in 1:1:360.
pred.sim(pred.sim==0)=360;
pred.mean(pred.mean==0)=360;

%draw
function drawRawPred5(pred,data,d,coh,pstd)
%data(cartesians)
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
            
%             %raw data
%             hold all
%             subplot(hs(1))
%             scatter(hs(1),d(F.inter.pos{k,i,j}),rawData{k,i,j}',...
%                 'MarkerEdgeColor',egdec(i,:),...
%                 'MarkerFaceColor',F.f2.color(i,:),...
%                 'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
% %             alpha(0.5)
% %             drawnow
%             
            %raw pred
            hold all
            subplot(hs(2))
            %colPred=F.f2.color(i,:)+0.2*([0 0 0]-F.f2.color(i,:));
            scatter(hs(2),d(F.inter.pos{k,i,j}),rawPred{k,i,j}',...
                'MarkerEdgeColor',egdec(i,:),...
                'MarkerFaceColor',F.f2.color(i,:),...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))));
%             alpha(0.5)
%             drawnow
            
%             %ideal
%             subplot(hs(1))
%             plot(F.f1.L,F.f1.L,'k.','Markersize',4)
%             plot(F.f1.L,225*ones(31,1),'b.','Markersize',4)
%             xlim([0 360])
%             ylim([0 360])
%             ylabel(hs(1),'Estimated directions (degrees)')
%             
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


%Von Mises
function mPdfs=vmPdfs(x,u,k,type)
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
%Calculate circular convolutions with circular data
function cconv=circConv(v1,v2)
cconv=ifft(fft(v1).*fft(v2));