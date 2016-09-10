
%usage
%[d upo]=simulateBasicData(0.06,0.1)

%How is a singe trial estimate generated in the brain?
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

function [d,upo]=simulateBasicData(km,kp,colori)

tic
%displayed direction. This should be replaced by the exact displayed
%direction we used in our task.
d=random('norm',225,80,5000,1);
kl=km;
up=225;

%make estimates
x=1:1:360;
upo=nan(numel(d),1);
for i=1:numel(d)
    
    %measurement disributions.
    m=vmPdfs(x,d(i),km);
   
    %sample
    mi=randsample(x,1,true,m);

    %trial-estimate.
    upo(i)=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)),(kl/kp)+cos(de2r(up,1)-de2r(mi,1))));    
end





%figure('color','w')
hold all
%simulation
plot(d,upo,'.','color',colori)
%ideal
plot(d,d,'--','color',[.7 .7 .7])
%prior's mean
plot(1:1:360,up(ones(360,1),:),'--b')

xlabel('Displayed directions(?)','fontsize',14)
ylabel('Simulated estimates(?)','fontsize',14)
xlim([0 361])
ylim([0 361])

toc






%NESTED
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