
% SLsimAttractor
% Pouget et al. 2011, PNAS


function SLsimAttractor

%%
numr=1000;

%response difference (rA-rB)
r=-numr/2:1:numr/2-1;
r=r/max(r);
%input (large gains for large r<0 and large r>0)
%it is zero when the two responses are equal
InputGain=[-numr/2:1:0 1:1:numr/2-1];
InputGain=InputGain/max(InputGain);
% InputGain=ones(1,numr)*100;

%white gaussian noise
wgn=randn(1,numr)*0.2;

%change of r with time
drdt=-4*r.*(r.^2 - 1) + InputGain + wgn;

%Energy landscape
Er=r.^2.*(r.^2 - 2) - InputGain.*r;

%draw
figure('color','w')
subplot(2,1,1)
axis square
hold all
plot(r,drdt,'r','linesmoothing','on')
plot([min(r) max(r)],[0 0],':','color',[.2 .2 .2])
ylabel('Change of r (drdt)')

subplot(2,1,2)
hold all
plot(r,Er,'r','linesmoothing','on')
plot([min(r) max(r)],[0 0],':','color',[.2 .2 .2])
xlabel('Response difference (r)')
ylabel('Energy')
axis square
SLremoveDeadSpace(0)