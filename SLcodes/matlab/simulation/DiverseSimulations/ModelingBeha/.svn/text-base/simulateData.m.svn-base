
%Usage:
    %k=[8.5 2.7 0.75 0.75 2.7 8.5 33];
    %n=10000;
    %[upo d coh pstd]=simulateData(1:10:360,k,n);

%Description:
    %k are the concentration parameters
    %n is the sample size

%Process: How is a singe trial estimate generated in the brain?
%Based on Girshick et al, 2010.
    %1.displayed direction.
    %2.measurement disribution.
    %mean is displayed direction, d, width is arbitrary (fit parameter)
    %3.sample a measurement, mi.
    %4.evaluate likelihood.
    %mean is the measurement mi, width is the width of the
    %measurement distribution kl=km;
    %5.prior.
    %experimental mean, width is arbitrary (fit parameter)
    %6.combine likelihood and prior and get trial-estimate (MAP).
    %mean of the posterior


function [upo,d,coh,pstd]=simulateData(k,n)

tic
%simulation - TO REMOVE(test)
%...............................................................
%displayed direction. This should be replaced by the exact displayed
%direction we used in our task.
x=5:10:355;

%simulate prior
pstd=randsample([80 40 20 10],n,'true');
prior(:,1)=(1./(2*sqrt(2*pi))).*exp(-((x-225).^2)./(2.*80^2));
prior(:,2)=(1/(2*sqrt(2*pi))).*exp(-((x-225).^2)/(2.*40^2));
prior(:,3)=(1/(2*sqrt(2*pi))).*exp(-((x-225).^2)/(2.*20^2));
prior(:,4)=(1/(2*sqrt(2*pi))).*exp(-((x-225).^2)/(2.*10^2));
Z=sum(prior);
prior=prior./Z(ones(size(prior,1),1),:);
pI(pstd==80)=1;
pI(pstd==40)=2;
pI(pstd==20)=3;
pI(pstd==10)=4;

%simulate directions
for i=1:numel(pstd)
    d(i)=randsample(x,1,'true',prior(:,pI(i)));
end

%simulate coherence
coh=randsample([0.24 0.12 0.06],numel(d),'true');
%...............................................................


%widths of the measurement distributions and priors
up=225;
kml(coh==0.24)=k(1);
kml(coh==0.12)=k(2);
kml(coh==0.06)=k(3);
kp(pstd==80)=k(4);
kp(pstd==40)=k(5);
kp(pstd==20)=k(6);
kp(pstd==10)=k(7);

%produce estimates
xe=1:1:360;
upo=nan(numel(d),1);

%measurement distribution
% m=vmPdfs(xe,d',kml);
% for i=1:numel(d)
%    
%     %sample
%     mi(i)=randsample(xe,1,true,m(:,i));
% end
% %trial-estimate.
% upo=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)), (kml./kp)+cos(de2r(up,1)-de2r(mi,1))));

%%trial-average prediction
% m=vmPdfs(xe,d',kml);
mi=d;
upo=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)), (kml./kp)+cos(de2r(up,1)-de2r(mi,1))));

%make sure data ranges in 1:1:360.
upo(upo==0)=360;

%PLOTS RAW
%figure('color','w')
figure('color','w')
subplot(1,3,1)
hold all
title('Coh:0.24')
%weak to strong
plot(d(pstd==80&coh==0.24),upo(pstd==80&coh==0.24),'.','color',[0.8 0.1 0])
plot(d(pstd==40&coh==0.24),upo(pstd==40&coh==0.24),'.','color',[1 0.4 0])
plot(d(pstd==20&coh==0.24),upo(pstd==20&coh==0.24),'.','color',[1 0.8 .2])
plot(d(pstd==10&coh==0.24),upo(pstd==10&coh==0.24),'.','color',[0.4 0.6 0])
plot(d,d,'--','color',[.7 .7 .7])
plot(1:1:360,up(ones(360,1),:),'--b')
xlabel('Displayed directions(?)')
ylabel('Estimated directions(?)')
xlim([0 361])
ylim([0 361])

subplot(1,3,2)
hold all
title('Coh:0.12')
plot(d(pstd==80&coh==0.12),upo(pstd==80&coh==0.12),'.','color',[0.8 0.1 0])
plot(d(pstd==40&coh==0.12),upo(pstd==40&coh==0.12),'.','color',[1 0.4 0])
plot(d(pstd==20&coh==0.12),upo(pstd==20&coh==0.12),'.','color',[1 0.8 .2])
plot(d(pstd==10&coh==0.12),upo(pstd==10&coh==0.12),'.','color',[0.4 0.6 0])
plot(d,d,'--','color',[.7 .7 .7])
plot(1:1:360,up(ones(360,1),:),'--b')
xlabel('Displayed directions(?)')
ylabel('Estimated directions(?)')
xlim([0 361])
ylim([0 361])

subplot(1,3,3)
hold all
title('Coh:0.06')
plot(d(pstd==80&coh==0.06),upo(pstd==80&coh==0.06),'.','color',[0.8 0.1 0])
plot(d(pstd==40&coh==0.06),upo(pstd==40&coh==0.06),'.','color',[1 0.4 0])
plot(d(pstd==20&coh==0.06),upo(pstd==20&coh==0.06),'.','color',[1 0.8 .2])
plot(d(pstd==10&coh==0.06),upo(pstd==10&coh==0.06),'.','color',[0.4 0.6 0])
plot(d,d,'--','color',[.7 .7 .7])
plot(1:1:360,up(ones(360,1),:),'--b')
xlabel('Displayed directions(?)')
ylabel('Estimated directions(?)')
xlim([0 361])
ylim([0 361])

%PLOT circular stats
[means,stds]=drawCircStat(upo,d,coh,pstd);


%draw circular statistics
function [means,stds]=drawCircStat(data,d,coh,pstd)
%Inputs: 
    %data
    %3 vectors for 3 factors. Rows are level values.

%data(cartesians)
data=polar2cartesian(data',1);

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
        
%Make mean & std
c=colormap;
F.f2.color={[0.5 0 0],...
    [1 0.2 0],...
    [1 0.6 0],...
    [0.75 0.75 0]};

%mean
figure('color','w')
for j=1:F.f2.n
    subplot(1,3,j)
    for k=1:F.f1.n
        for i=1:F.f3.n
            stat{k,i,j}=statcircular(data(F.inter.pos{k,i,j},:));
            means(k,i,j)=stat{k,i,j}.deg.mean;
        
            %draw
            hold all
            scatter(F.f1.L(k),means(k,i,j)',...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))
            %ideal
            plot(F.f1.L(k),F.f1.L(k),'k.','Markersize',3)
        end
    end
    xlim([0 360])
    ylim([0 360])
end

%std
figure('color','w')
for j=1:F.f2.n
    subplot(1,3,j)
    for k=1:F.f1.n
        for i=1:F.f3.n
            stds(k,i,j)=stat{k,i,j}.deg.std;       
        
            %draw
            hold all
            scatter(F.f1.L(k),stds(k,i,j)',...
                'MarkerEdgeColor','w',...
                'MarkerFaceColor',F.f2.color{i},...
                'displayname',strcat(F.f2.nm,':',num2str(F.f2.L(j))))

        end
    end
    xlim([0 360])
    ylim([0 360])
end



toc
%von Mises pdf
function mPdfs=vmPdfs(x,u,k)
%k must be in the range [0:709]. when k>709, probabilities are inf.
%k must be in the range [0:700]. when k>700, besseli is inf amd mPdfs is 0.
%When normalizing the von mises to probabilities we get NaN values
%(0/sum(0)).
w=1;

%radians
x=de2r(x,1); x=x';
u=de2r(u,1); u=u';
x2=x(:,ones(numel(u),1));
u2=u(ones(numel(x),1),:);
k2=k(ones(numel(x),1),:);
mPdfs=exp(k2.*cos(w*(x2-u2)))./(2*pi.*besseli(0,k2));
%pdf
Z_=sum(mPdfs);
Z=Z_(ones(numel(x),1),:);
mPdfs=mPdfs./Z;
%degrees to radians
function radians=de2r(ang,sign)

%not signed radians (1:2*pi)
radians=(ang/360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end
%convert radians to degrees
function degrees=ra2d(theta)

%When input radians are between 0:2*pi
degrees=(theta/(2*pi))*360;

%if larger than 360 degrees then subtract
%360 degrees
while (sum(degrees>360))
    degrees = degrees - (degrees>360)*360;
end

%if less than 360 degreees then add
%360 degrees
while (sum(degrees<-360))
    degrees = degrees + (degrees<-360)*360;
end

%When radians are signed between -pi:pi.
degrees(degrees<0)=degrees(degrees<0)+360;
%Calculate the circular statistics of the data
function data=statcircular(coord)
%input a vector of cartesian coordinates (coord)

%register the coordinates of the input directions
data.coord.all=coord;

%convert from cartesian coordinates to angles (in degree)
data.deg.all=getangle(coord(:,1),coord(:,2));

%calculate the cartesian coordinates of the mean direction est
data.coord.mean=nanmean(coord,1);

%calculate the mean direction estimate (in degree)
data.deg.mean=getangle(data.coord.mean(:,1),data.coord.mean(:,2));

%calculate the std to the mean direction est (in degree);!!! could be a
%subfunction itself
%Apply the rule of thumb that follows. It seems to work fine intuitively. It would be nice to
%fine a cleaner way to calculate the std.
%initialize the 'sample' and 'mean' variables used to calculate the std
data.num=numel(data.deg.all);%sample size
data.deg.allforstd=data.deg.all;
data.deg.meanforstd=repmat(data.deg.mean,data.num,1);

%if the resulting mean direction is between 0 and 180.
if data.deg.mean+180<=360
    %if estimation is>=mean direction+180
    for i=1:data.num
        if data.deg.all(i)>=data.deg.mean+180
            data.deg.allforstd(i)=data.deg.all(i)-360;
        end
    end
    %if the resulting mean direction is between 180 and 360.
else
    %if the estimated direction sampled is <=the mean direction - 180
    for i=1:data.num
        if data.deg.all(i)<=data.deg.mean-180
            data.deg.meanforstd(i)=data.deg.mean-360;
        end
    end
end

%now calculate the variance of the estimated direction.
data.deg.var=nanmean((data.deg.allforstd-data.deg.meanforstd).^2,1);

%and now calculate the std
data.deg.std=sqrt(data.deg.var);

%and now calculate the sem
data.deg.sem=data.deg.std/sqrt(data.num);
%polar to cartesian
function v=polar2cartesian(di,r)
%di:angs in degree
%r:radius of the unit circle

%degrees & radians
di2.deg=di;
di2.rad=(di2.deg*pi)/180;
% di2.rad=(di2.deg*3.14)/180;


%cartesian vectors
x=r.*cos(di2.rad);
y=r.*sin(di2.rad);
v=[x y];
%Convert from cartesian coordinates to angles (degree)
function [angle]=getangle(x,y)
% check! to check if the function works fine, write
% e.g., input=180; output=getangle(cos(angle*pi/180),sin(angle*pi/180));
% if the function works input and output should always be the same between
% 0 and 360.

% convert from cartesian coordinates to angle in radians
angle = atan(y./x); % theta = arctan(opposite/adjacent);

% adjust each angle according to his quadrant (in degree)
for i=1:numel(x) % sample each angle
    if x(i)>=0 && y(i)>=0                   %(quadrant 1)
        angle(i) = angle(i)*180/pi;
    elseif x(i)<0                      %(quadrant 2 & 3)
        angle(i) = angle(i)*180/pi + 180;
    elseif x(i)>=0 && y(i)<0               %(quadrant 4)
        angle(i) = angle(i)*180/pi + 360;
    end
end