
function [upo]=MakePre5(d,coh,pstd,k)
tic
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

%%measurement distribution
%m=vmPdfs(xe,d',kml);
%for i=1:numel(d)
%    
%     %sample
%     mi(i)=d(i);%randsample(xe,1,true,m);
%end
%%trial-estimate.
%upo=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)), (kml./kp)+cos(de2r(up,1)-de2r(mi,1))));

%trial-average prediction
m=vmPdfs(xe,d',kml);
mi=d;
upo=ra2d(de2r(mi,1)+atan2(sin(de2r(up,1)-de2r(mi,1)),(kml./kp)+cos(de2r(up,1)-de2r(mi,1))));
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
mPdfs=exp(k2.*cos(w.*(x2-u2)))./(2*pi.*besseli(0,k2));
%pdf
Z_=sum(mPdfs);
Z=Z_(ones(numel(x),1),:);
mPdfs=mPdfs./Z;
%degrees to radians
function radians=de2r(ang,sign)

%not signed radians (1:2*pi)
radians=(ang./360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end
%convert radians to degrees
function degrees=ra2d(theta)

%When input radians are between 0:2*pi
degrees=(theta./(2*pi))*360;

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
