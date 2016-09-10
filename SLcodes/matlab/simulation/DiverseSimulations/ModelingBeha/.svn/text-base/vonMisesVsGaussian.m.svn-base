%vonMisesVsGaussian
figure('color','w'); 
hold all; 
vm=vmPdfs(1:1:360,225,2.7714,'norm');
f=(1/(s*sqrt(2*pi)))*exp(-((x-mu).^2)/(2*s.^2));
f=f/sum(f);
plot(vm,'r','linesmoothing','on'); 
plot(f2,'k','linesmoothing','on')
set(gca,'fontsize',14)
xlabel('Displayed direction'); ylabel('Probability')
legend('Von mises (k=2.7714)','Gaussian(std=40)')