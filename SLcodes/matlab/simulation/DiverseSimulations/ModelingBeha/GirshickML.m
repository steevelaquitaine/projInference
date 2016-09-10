
%usage: input data as upo.
    %e.g., [logL ML]=GirshickML(0,225,10,1)

%Description
    %Inputs can be vectors where rows are associated with each others.
    %e.g, data 0 is explained by displayed direction 225, a kpi=10 and
    %kli=1 with likelihood logL.
    
function [logL,ML]=GirshickML(data,displ,kpi,kli)
%Calculate the likelihood of observing any given data with Girshick estimator
%model.
%kpi and kli cannot both be 0.(kli/kpi) will collapse to NaN in the bayesian
%inference equation that gives the mean of the posterior.

%parameters (degrees)
upo=data;
up=225; 
ulall=1:1:360;

%method1: ul from Bayes closed form equation (vm)
%options=optimset('TolFun',10^-10);
%ul=fsolve(@(ul) upo-(ul+atan2(sin(up-ul),(kli/kpi)+cos(up-ul))),pi,options);

%method2:row upo/col ul; fun is the function to minimize to find ul
%solution(degrees).ulall2, up and upo must be in signed radians because the von
%mises from which the equation is derived takes radians values for ulall2,
%and up. 
upo2=upo(:,ones(1,numel(ulall)));
ulall2=ulall(ones(numel(upo),1),:);
ki=(kli./kpi)';
ki=ki(:,ones(1,numel(ulall)));
fun=abs(de2r(upo2,1)-(de2r(ulall2,1)+atan2(sin(de2r(up,1)-de2r(ulall2,1)),ki+cos(de2r(up,1)-de2r(ulall2,1)))));
[dummy,I]=min(fun,[],2);
%degrees
ul=ulall(I);

%prior
x=ulall;%1:1:360;
up=up*ones(numel(upo),1);
p=vmPdfs(x,up,kpi);
Z=sum(p);
p=p./Z(ones(numel(x),1),:);


%posterior(approximation)
kpoi=sqrt((kli.^2)+(kpi.^2)+(2*kli.*kpi.*cos(de2r(up,1)'-de2r(ul,1))));
po.apx=vmPdfs(x,upo,kpoi);
Z=sum(po.apx);
po.normapx=po.apx./Z(ones(numel(x),1),:);
[~,I]=max(po.normapx);
MAP.normapx=x(I);

%posterior(exact)
W=(besseli(0,kpoi)./(2*pi*besseli(0,kli).*besseli(0,kpi)));
po.cf=po.apx.*W(ones(numel(x),1),:);
Z=sum(po.cf);
po.cf=po.cf./Z(ones(numel(x),1),:);
[~,I]=max(po.cf);
MAP.cf=x(I);

%llh (what we tried to solved for)
l=vmPdfs(x,ul',kli);
Z=sum(l);
l=l./Z(ones(numel(x),1),:);

%The likelihood of observing data i is the probability that measurement
%mi, that is the mean of the likelihood (ul(i)), is observed given a
%measurement distribution which mean is the displayed direction.
um=displ;
m=ulall;%1:1:360;
mpdf=vmPdfs(m,um',kli);
Z=sum(mpdf);
mpdf=mpdf./Z(ones(numel(x),1),:);

%Single trial's measurement, its position(row) for each trial(col) and its
%probability (also likelihood of trial's data). The probability of
%observing a particular "mi" given the displayed direction "displ" is equal 
%to the maximum of the likelihood function that peaks at mi. Note that 
%likelihood is not normalized. Likelihood is combined with the prior to 
%give a posterior density. The mode of the posterior (MAP estimate) is the 
%estimated direction.
mi=ul;
%vectorized way
%mipos=mi'+numel(m).*((0:1:numel(mi)-1)');
%ML.p=mpdf(mipos);
%loop way
for i=1:numel(mi)
    ML.p(i)=mpdf(ulall==mi(i),i);
end

%ul-true
%7.4609
%mi-true
%0.1302
%ul-fit
%10.55
%mi-fit
%0.18



%adjusted likelihood.
%lapserate=1e-10;
%ML.p=ML.p+lapserate;

%adjusted log likelihood
logL=log(ML.p);

%-SumLogL
minsumlogL=-sum(logL);

% %plot
% figure('color','w')
% hold all
% %measurement density
% plot(x,mpdf,'r')
% %likelihood of measurement i origin of the observed MAP
% plot(x,l','k')
% %prior
% plot(x,p','--k')
% %posterior
% plot(x,po.normapx','b')
% plot(x,po.cf',':y')
% %measurement i
% plot([ML.mi ML.mi],[0 ML.p],'.r','markersize',20)
% %MAP
% plot([MAP.normapx MAP.normapx],[0 max(po.normapx)],'--b')
% %plot([MAP.cf MAP.cf],[0 max(po.cf)],':y','linewidth',1.5)
% %annotations
% text(ML.mi,ML.p+0.001,0,num2str(ML.mi),'color','r')
% text(ul,max(l),0,num2str(ul),'color','k')
% text(up,max(p),0,num2str(up))
% text(upo,max(po.normapx),0,num2str(upo),'color','b');
% text(um,max(mpdf)/2,0,'measurement pdf','color','r')
% text(ul,max(l)/2,0,'likelihood','color','k')
% text(up,max(p)/2+0.01,0,'prior')
% text(upo,max(po.normapx)/2,0,'posterior','color','b');
% xlim([0 360])




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
%Z_=sum(mPdfs);
%Z=Z_(ones(numel(x),1),:);
%mPdfs=mPdfs./Z;
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

