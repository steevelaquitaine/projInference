
%steeve laquitaine

%Usage
    %[S1 S2]=simulateVMStd(1:1:360,180,300);

% Stdev1: std of simulated samples
% Stdev2: std of discrete von mises distribution
        
%The maximum standard deviation measured in my data was 100 degrees for
%directions most distant to prior's mean. This corresponds to a k of 0
%degrees corresponding to a flat von Mises. 

%The sharpest prior condition produced data with std~10 degrees. This 
%corresponds to k~33.

%k=699 corresponds to a std of about 2 degrees.

function [Stdev1 Stdev2]=simulateVMStd(x,u,k)
%VMSTD Summary of this function goes here
%Detailed explanation goes here
n=1000;
p=vmPdfs(x,u,k,'norm');
y.deg=randsample(x,n,true,p);


%histogram samples from von Mises.
psampled=hist(y.deg,x)./sum(hist(y.deg,x));
figure('color','w');
hold all
bar(x,psampled,'w','edgecolor',[.7 .7 .7]);
plot(p,'r','linewidth',2)
xlim([0 360])

%two ways of calculating std
Stdev1=std(y.deg);
Stdev2=sqrt(sum(p.*((x-u).^2)'));

%Draw corresponding gaussian
plot(normpdf(x,u,Stdev1),':','color',[0 0.5 1],'linewidth',2)
xlabel('x')
ylabel('Probability')
ylim([0 0.17])
box off


%Draw circular
%convert degrees in polar angle (rad,radius)
figure('color','w')
y.rad=de2r(y.deg,0);
rho=ones(n,1)*2.5;
polar(y.rad',rho,'o');

%von mises
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
%degrees to radians
function radians=de2r(ang,sign)

%not signed radians (1:2*pi)
radians=(ang/360)*2*pi;

%sign radians(-pi:pi)
if sign==1
    radians(ang>180)=(ang(ang>180)-360)*(2*pi/360);
end
