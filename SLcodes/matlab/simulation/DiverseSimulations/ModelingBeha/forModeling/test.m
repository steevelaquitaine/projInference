
%This setting permits to obtain a good approximation of the continous
%posterior obtained from bayesian inference when multiplying two discretes
%pdfs. The "sample" (Jazayeri) space x must be large enough so that the
%discrete pdfs approximate properly their continuous version.

clear; close all

uLl=[5:10:355]';
uPr=225;
%x must be large enough.
x=[-1000:10:1000]';
%sl and sp must be small enough. The approximation is exact below 100.
sl=ones(36,1)*80;
sp=ones(36,1)*80;
sM=0.01;

for i=1:numel(uLl)
    %llh
    l(:,i)=((sl(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uLl(i)).^2)/(sl(i).^2));
    l(:,i)=l(:,i)/sum(l(:,i));    
    
    %prior
    pr(:,i)=((sp(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uPr).^2)/(sp(i).^2));
    pr(:,i)=pr(:,i)/sum(pr(:,i));
    
    %posterior
    po(:,i)=l(:,i).*pr(:,i)/sum(l(:,i).*pr(:,i));
    uPo(i,1)=sum(x.*po(:,i));
    sPo(i,1)=sum(sqrt(po(:,i).*(x-uPo(i)).^2));
    
    %estimate
    sEs(i,1)=sPo(i)+sM;
   
    hold all
    plot(x,l(:,i))
    plot([uLl(i) uLl(i)],[0 0.01],'--k')
    plot([sum(x.*l(:,i)) sum(x.*l(:,i))],[0 0.01],'--r')

    drawnow
    pause(0.)
end


for i=1:numel(uLl)
    ucheck=sum(x.*l(:,i));
    uLl(i);
    fprintf('\n %f %f \n',ucheck,uLl(i))
end


%% How to make the fitting quicker?
clear; close all

uLl=[5:10:355]';
uPr=225;
%sl and sp must be small enough. The approximation is exact below 100.
sl=ones(36,1)*80;
sp=ones(36,1)*80;
sM=0.01;






%%
l=@(x) ((sl*sqrt(2*pi)).^-1).*exp(-0.5*((x-uLl).^2)/(sl.^2));
pr=@(x) ((sp*sqrt(2*pi)).^-1).*exp(-0.5*((x-uPr).^2)/(sp.^2))/sum(((sp*sqrt(2*pi)).^-1).*exp(-0.5*((x-uPr).^2)/(sp.^2)));
upo=@(x) x.*(l(x).*pr(x))/sum(l(x).*pr(x));
integral(upo,-1000,1000);


integral(l,-1000,1000)

%%





for i=1:numel(uLl)
    %llh
    l(:,i)=((sl(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uLl(i)).^2)/(sl(i).^2));
    l(:,i)=l(:,i)/sum(l(:,i));    
    
    %prior
    pr(:,i)=((sp(i)*sqrt(2*pi)).^-1).*exp(-0.5*((x-uPr).^2)/(sp(i).^2));
    pr(:,i)=pr(:,i)/sum(pr(:,i));
    
    %posterior
    po(:,i)=l(:,i).*pr(:,i)/sum(l(:,i).*pr(:,i));
    uPo(i,1)=sum(x.*po(:,i));
    sPo(i,1)=sum(sqrt(po(:,i).*(x-uPo(i)).^2));
    
    %estimate
    sEs(i,1)=sPo(i)+sM;
   
    hold all
    plot(x,l(:,i))
    plot([uLl(i) uLl(i)],[0 0.01],'--k')
    plot([sum(x.*l(:,i)) sum(x.*l(:,i))],[0 0.01],'--r')

    drawnow
    pause(0.)
end


for i=1:numel(uLl)
    ucheck=sum(x.*l(:,i));
    uLl(i);
    fprintf('\n %f %f \n',ucheck,uLl(i))
end
