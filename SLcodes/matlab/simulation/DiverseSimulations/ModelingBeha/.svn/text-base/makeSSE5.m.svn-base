
%Alternative to GirshickML

function [logL,fitP]=makeSSE5(data,displ,coh,Pstd,fitP)

%This nethod is 0.20 seconds faster per function evaluation than 1)
%re-calculating the solution for each data point in a loop and 2) >20
%seconds faster per function evaluation than using fsolve to find the
%solution for each data point in the loop.

%This code can be improved look at GirshickML code in Modeling folder. We
%can get rid of the rounding of the data and of the for loop!


%fit parameters(k has no unit, ideally k is in the range 0:inf)
kl1=fitP(1);
kl2=fitP(2);
kl3=fitP(3);
kp1=fitP(4);
kp2=fitP(5);
kp3=fitP(6);
kp4=fitP(7);

%degrees
ulall=1:1:360;
up=225;
ul=nan(numel(data),1);

%3coh*4priors Girshick lookup tables (degrees)
ul11=lookuptable(ulall,up,kl1,kp1);
ul12=lookuptable(ulall,up,kl1,kp2);
ul13=lookuptable(ulall,up,kl1,kp3);
ul14=lookuptable(ulall,up,kl1,kp4);
ul21=lookuptable(ulall,up,kl2,kp1);
ul22=lookuptable(ulall,up,kl2,kp2);
ul23=lookuptable(ulall,up,kl2,kp3);
ul24=lookuptable(ulall,up,kl2,kp4);
ul31=lookuptable(ulall,up,kl3,kp1);
ul32=lookuptable(ulall,up,kl3,kp2);
ul33=lookuptable(ulall,up,kl3,kp3);
ul34=lookuptable(ulall,up,kl3,kp4);

%assign ul solution (degrees) to each data point
data=round(data)';
%parpool parfor if possible
%make sure data ranges in 1:1:360.I saw that simulateData produces both 0
%and 360 that are the same. All 0 will be 360.
data(data==0)=360;
for i=1:360
    ul(data==i&coh==0.24&Pstd==80)=ul11(i);
    ul(data==i&coh==0.24&Pstd==40)=ul12(i);
    ul(data==i&coh==0.24&Pstd==20)=ul13(i);
    ul(data==i&coh==0.24&Pstd==10)=ul14(i);
    ul(data==i&coh==0.12&Pstd==80)=ul21(i);
    ul(data==i&coh==0.12&Pstd==40)=ul22(i);
    ul(data==i&coh==0.12&Pstd==20)=ul23(i);
    ul(data==i&coh==0.12&Pstd==10)=ul24(i);
    ul(data==i&coh==0.06&Pstd==80)=ul31(i);
    ul(data==i&coh==0.06&Pstd==40)=ul32(i);
    ul(data==i&coh==0.06&Pstd==20)=ul33(i);
    ul(data==i&coh==0.06&Pstd==10)=ul34(i);
end

%The likelihood of observing data i is the probability that measurement
%mi, that is the mean of the likelihood (ul(i)), is observed given a
%measurement distribution which mean is the displayed direction.
%degrees
m=1:1:360;
um=displ';
kl(coh==0.24)=fitP(1);
kl(coh==0.12)=fitP(2);
kl(coh==0.06)=fitP(3);
mpdf=vmPdfs(m,um,kl);

%Single trial's measurement, its position(row) for each trial(col) and its
%probability (also likelihood of trial's data). The probability of
%observing a particular "mi" given the displayed direction "displ" is equal 
%to the maximum of the likelihood function that peaks at mi. Note that 
%likelihood is not normalized. Likelihood is combined with the prior to 
%give a posterior density. The mode of the posterior (MAP estimate) is the 
%estimated direction.
mi=ul;
% mipos=mi+numel(m)*((0:1:numel(mi)-1)');
% ML.p=mpdf(mipos);
for i=1:numel(mi)
    ML.p(i)=mpdf(ulall==mi(i),i);
end
%adjusted log likelihood.
%lapserate=1e-10;
%ML.p=ML.p+lapserate;

%We use log likelihood because likelihood is so small that matlab cannot
%encode it properly (numerical unstability).
logL=-sum(log(ML.p));
%logL=-sum(ML.p);

%Girshick lookup table
function ul=lookuptable(ulall,up,kl,kp)
%method2:row upo/col ul; fun is the function to minimize to find ul
%solution(degrees).ulall2, up and upo must be in signed radians because the von
%mises from which the equation is derived takes radians values for ulall2,
%and up.
%lookuptable inputs are in degrees and are converted in radians
% upo=[1:0.1:360]';
upo=[1:1:360]';

upo2=upo(:,ones(1,numel(ulall)));
ulall2=ulall(ones(numel(upo),1),:);
fun=abs(de2r(upo2,1)-(de2r(ulall2,1)+atan2(sin(de2r(up,1)-de2r(ulall2,1)),(kl/kp)+cos(de2r(up,1)-de2r(ulall2,1)))));
[dummy,I]=min(fun,[],2);

%in degrees
ul=ulall(I);