

%SNR=stdFR; SNR=stdllh; stdllh=stdFR;
%stdllh=sqrt(1/c);stdFR=sqrt(1/c);varFR=1/c;
clf
c=0:0.001:1;
varllh=1./c;
stdllh=sqrt(1./c);

%plot stdllh vs. c.
subplot(1,2,1)
plot(c,stdllh,'-','linewidth',2)
xlim([0 1])
ylim([0 max(1./sqrt(c))])
xlabel('coherence')
ylabel('std llh')
box off

%plot varFR vs. c.
subplot(1,2,2)
plot(c,1./c,'-','linewidth',2)
xlim([0 1])
ylim([0 max(1./sqrt(c))])
xlabel('coherence')
ylabel('variance FR')
box off

figure;
plot(c,1.1.*c)