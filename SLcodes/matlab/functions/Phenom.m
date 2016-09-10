
% steeve laquitaine

% need a vector of data : P=basic fitting polynom 10 (fx)
% give the derivative function 2 of t+1 - t






P=fx;

for t=1:length(P)-1  
    dP(t+1)=P(t+1)-P(t)
end

dP=dP.*dP; % dP au carré meilleur visualisation des derivé diff de 0
dP=dP';

x=1:length(P);
tresh=ones(length(P),1)*0.0000025  %0.05 au carré


% dynamic phase
i1=find(dP>=0.0000025)
LP(1:length(dP),1)=NaN;
LP(i1)=dP(i1);            
phs=max(i1)



plot(dP,'.');
hold on; plot(x,tresh,'r')
hold on; plot(LP,'g')