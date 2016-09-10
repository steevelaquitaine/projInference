
close all

subplot(121)
llh = vmPdfs(1:1:360,135,20,'norm');
prior = vmPdfs(1:1:360,225,0,'norm');
pos = llh.*prior/sum(llh.*prior);

hold all
plot(llh)
plot(prior)
plot(pos,'r')
ylim([0 0.05])


subplot(122)
llh = vmPdfs(1:1:360,135,20,'norm');
prior = vmPdfs(1:1:360,225,20,'norm');
pos = llh.*prior/sum(llh.*prior);

hold all
plot(llh)
plot(prior)
plot(pos,'r')
ylim([0 0.05])
