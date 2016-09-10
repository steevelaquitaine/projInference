
 %author: steeve laquitaine
%purpose: get the probability of prior switch based on softmax comparison
%of prior and llh strengths.

%Effect of Beta for weak coherence and prior.
%The range of prior strengths covers the experimental priors. The llh
%strength is intermediate

function [Ppriormn,Pllhmn,Prnew]=simHeuristicSoftmax(kp,kml,Pr,Beta)

%softmax competition
for i=1:numel(Beta)
    weightPriormn(:,i)=exp(Beta(i).*log(kp))./(exp(Beta(i).*log(kp))+exp(Beta(i).*log(kml)));
end
weightLlhmn=1-weightPriormn;

%scale "weights" to probabilities
sumP=weightPriormn+weightLlhmn+Pr;
Ppriormn=1-(Pr+weightLlhmn)./sumP;
Pllhmn=1-(Pr+weightPriormn)./sumP;
Prnew=1-(weightLlhmn+weightPriormn)./sumP;

%draw
figure('color','w')
plot(Ppriormn,'linewidth',1.001,'linesmoothing','on')
box off
set(gca,'fontsize',14)
xlabel({'Prior strength (kp) increases','llh strength=2'})
ylabel('Probability of prior switch')
xlim([kp(1) kp(end)])