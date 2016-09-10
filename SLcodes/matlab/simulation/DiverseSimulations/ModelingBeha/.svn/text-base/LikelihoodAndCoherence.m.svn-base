

%The width of the likelihood depends on the SNR, MeanFR/stdFR.
%SNR=mean/sqrt(var)
%mean=var; 
%--> so SNR=(mean*sqrt(var))/var=sqrt(var);
%--> SNR=std

%If SNR is the width km of the measurement distribution, 
%p(m/di), for a motion direction (ui). It is also the 
%width of the likelihood, p(mi,d), of 
%observing any measurement mi sampled from this measurement
%distribution.

%Variance FR linearly increases with coherence:
%var=coh. So SNR=sqrt(coh) meaning that km decreases linearly with sqrt
%(coh).

%britten et al, 1993
b=1.2;
k=1;
varFR=k.*meanFR.^b;

%SNR increases linearly, varllh decreases linearly.
%meanFR/sqrt(varFR) increases linearly, varllh decreases.
%coh/sqrt(coh), i.e., sqrt(coh) increases, varllh decreases.
hold all
coh=0:0.1:1;
plot(coh,sqrt(coh))
varllh=0.1./coh;
plot(coh,varllh,'r')



