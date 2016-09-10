
    %Author: steeve Laquitaine
      %date: 131024
      %usage:
        %load data03
        %k=[80 40 15 0.74559 2.77 8.74 33.25 10.^-10000 15 0.07];%sub01 best
        %pred=drawPred03(data,disp,coh,pstd,k);

%Note: One can use the code with "drawnow" but it's 6 times faster without.

%Note: careful: I couldn't implement this code.
%a single "upo" must always correspond to a single "mi(ul)"
%according to bayesian inference. So we should not find duplicate "mi".
%It seems that it is not possible to start from each possible upo and
%reverse back to its corresponding "mi" such that we get unique values of 
%"mi". We would need to increase the resolution of ulall which would
%significantly slow down the code. I also noticed that for high kp (10000)
%the equation in the lookup table produces aberrant "ul". It is difficult
%to manipulate systematically the range of values that don't produces
%aberrant values.
%Need to think about another approach

function pred=drawPred03(data,d,coh,pstd,k)
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
xe=1:1:360;


%P(estimate/BI)
%---------------
%parameters (degrees)
kl1=k(1);
kl2=k(2);
kl3=k(3);
kp1=k(4);
kp2=k(5);
kp3=k(6);
kp4=k(7);
ulall=1:1:360;
up=225;
upo=[1:1:360]';
ulall=1:1:360;

%Find the "ul" that produced "upo" ranging from 1:360 for the 12 
%conditions (3 X 4 lookup tables(degrees))
%ul11(1) produced upo(1) and ul11(2) produced upo(2), etc...
ul11=lookuptable(upo,ulall,up,kl1,kp1);
ul12=lookuptable(upo,ulall,up,kl1,kp2);
ul13=lookuptable(upo,ulall,up,kl1,kp3);
ul14=lookuptable(upo,ulall,up,kl1,kp4);
ul21=lookuptable(upo,ulall,up,kl2,kp1);
ul22=lookuptable(upo,ulall,up,kl2,kp2);
ul23=lookuptable(upo,ulall,up,kl2,kp3);
ul24=lookuptable(upo,ulall,up,kl2,kp4);
ul31=lookuptable(upo,ulall,up,kl3,kp1);
ul32=lookuptable(upo,ulall,up,kl3,kp2);
ul33=lookuptable(upo,ulall,up,kl3,kp3);
ul34=lookuptable(upo,ulall,up,kl3,kp4);

%calculate p(estimate/BI,d,km,kp) in each condition
%First, calculate measurement densities produced by motion direction at
%each trial.
%Note: careful: a single "upo" must always correspond to a single "mi(ul)"
%according to bayesian inference. So we should not find duplicate "mi".
%It seems that it is not possible to start from each possible upo and
%reverse back to its corresponding "mi" such that we get unique values of 
%"mi". We would need to increase the resolution of ulall which would
%significantly slow down the code. I also noticed that for high kp (10000)
%the equation in the lookup table produces aberrant "ul". It is difficult
%to manipulate systematically the range of values that don't produces
%aberrant values.
%Need to think about another approach
p_est=nan(360,numel(d));
p_meas=vmPdfs(xe,d,kml,'norm');
for triali=1:numel(d)
    if kml(triali)==kl1 && kp(triali)==kp1
        mi=ul11;
    elseif kml(triali)==kl1 && kp(triali)==kp2
        mi=ul12;
    elseif kml(triali)==kl1 && kp(triali)==kp3
        mi=ul13;
    elseif kml(triali)==kl1 && kp(triali)==kp4
        mi=ul14;
    elseif kml(triali)==kl2 && kp(triali)==kp1
        mi=ul21;
    elseif kml(triali)==kl2 && kp(triali)==kp2
        mi=ul22;
    elseif kml(triali)==kl2 && kp(triali)==kp3
        mi=ul23;
    elseif kml(triali)==kl2 && kp(triali)==kp4
        mi=ul24;
    elseif kml(triali)==kl3 && kp(triali)==kp1
        mi=ul31;
    elseif kml(triali)==kl3 && kp(triali)==kp2
        mi=ul32;
    elseif kml(triali)==kl3 && kp(triali)==kp3
        mi=ul33;
    elseif kml(triali)==kl3 && kp(triali)==kp4
        mi=ul34;
    end   
    for est=1:360
        p_est(est,triali)=p_meas(xe==mi(est),triali);
    end
end
toc


% KEEP CHECKING P_EST !!!!!!


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

preallocate
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

%Girshick lookup table
function ul=lookuptable(upo,ulall,up,kl,kp)
%fun is the function to minimize to find ul
%solution(degrees).ulall2, up and upo must be in signed radians because the von
%mises from which the equation is derived takes radians values for ulall2,
%and up. lookuptable inputs are in degrees and are converted in radians
%upo=1;kl=kl3;kp=kp1;%debug
upo2=upo(:,ones(1,numel(ulall)));
ulall2=ulall(ones(numel(upo),1),:);

%The "ulall2" that minimizes the function is the solution. One "upo" should 
%have a unique "ulall2" solution for a given pair [kl,kp].
%"fun"=0 when "de2r(upo2,1)" = right side of the equation and we get the
%solution to "ulall2".
%! note: careful because in some situations, e.g., with extremely large kp
%the equation produces aberrant "ul".
fun=abs(de2r(upo2,1) - (de2r(ulall2,1)+atan2(sin(de2r(up,1)-de2r(ulall2,1)),(kl/kp)+cos(de2r(up,1)-de2r(ulall2,1)))));
[dummy,I]=min(fun,[],2);
ul=ulall(I);

%degrees to radians
function radians=de2r(ang,sign)

%not signed radians (1:2*pi)
radians=(ang/360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end
