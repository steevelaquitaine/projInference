function [MAP ML]=BayesVM
%VMPRODUCT Summary of this function goes here
%Detailed explanation goes here
%The posterior that results from the product of two von mises densities is 
%not a von Mises (see Murray and Morgenstern, 2010, journal of vision).
%We check that the posterior resuting from the product of two von mises is
%well approximated by a von Mises. 
%The maximum likelihood estimates and probability to observe those estimates 
%were identical between the three cases.


%Simulation of Bayes' exact posterior
subplot(3,1,1)
x=[0:1:359];
ul=90;
up=0;
kl=42;
kp=42;

%llh
l=vmPdfs(x,ul,kl);
l=l/sum(l);

%prior
p=vmPdfs(x,up,kp);
p=p/sum(p);

%posterior
po.si=(l.*p)/sum(l.*p);

%MAP
[ML.si,I]=max(po.si);
MAP.si=x(I);

%plot
hold all
plot(x,po.si','k')
plot([MAP.si MAP.si],[0 max(po.si)],'--b')




%Closed form equation of Bayes' exact posterior
subplot(3,1,2)
ul=de2r(ul);
up=de2r(up);
upo=ul+atan2(sin(up-ul),(kl/kp)+cos(up-ul));
kpo=sqrt((kl^2)+(kp^2)+(2*kl*kp*cos(up-ul)));
upo=ra2d(upo);
vm.approx=vmPdfs(x,upo,kpo);
po.cf=vm.approx*(besseli(0,kpo)/(2*pi*besseli(0,kl)*besseli(0,kp)));
po.cf=po.cf/sum(po.cf);

%MAP
[ML.cf,I]=max(po.cf);
MAP.cf=x(I);

hold all
plot(x,po.cf,'color',[.5 .5 .5])
plot([MAP.cf MAP.cf],[0 max(po.cf)],'--r')





%Bayes' approximated posterior
subplot(3,1,3)
po.approx=vm.approx;
po.approx=po.approx/sum(po.approx);

%MAP
[ML.approx,I]=max(po.approx);
MAP.approx=x(I);

hold all
plot(x,po.approx,'color',[.5 .5 .5])
plot([MAP.approx MAP.approx],[0 max(po.approx)],'--r')







function mPdfs=vmPdfs(x,u,k)
w=1;
%radian
u=de2r(u); u=u';
x=de2r(x); x=x';
u2=u(ones(numel(x),1),:);
x2=x(:,ones(numel(u),1));
%mPdfs=nan(numel(x),numel(u));
mPdfs=exp(k.*cos(w*(x2-u2)))/2*pi.*besseli(0,k);
% %pdf
% Z_=sum(mPdfs);
% Z=Z_(ones(numel(x),1),:);
% mPdfs=mPdfs./Z;

%degrees to radians
function radians=de2r(ang)
radians=(ang/360)*2*pi;
%convert radians to degrees
function degrees=ra2d(theta)
degrees=(theta/(2*pi))*360;

% if larger than 360 degrees then subtract
% 360 degrees
while (sum(degrees>360))
  degrees = degrees - (degrees>360)*360;
end

% if less than 360 degreees then add 
% 360 degrees
while (sum(degrees<-360))
  degrees = degrees + (degrees<-360)*360;
end