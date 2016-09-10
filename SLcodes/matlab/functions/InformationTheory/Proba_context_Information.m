
% Mean Quantity of information (predictibility of event) by receveived
% reinforcement in the different uncertainty context


% Information Theory or uncertainty level: Shannon:
% H(x)=-sum(p(i)*log2(p(i))) with message x composed of 2 symbol: i=1 or
% i=2 and H(x)=H(x1)+ H(x2)

% i=1:reward;  i=0:noreward

%Infx1x2 = H(x1) + H(x2)= (-p(1)*log2(p(1) - p(0)*log2(p(0))) + (-p(0)*log2(p(0))-p(1)*log2(p(1))

% E(r1)*log2E(r1) – (1-E(r1))*log2E(1-E(r1)) - E(r2)*log2E(r2) – (1-E(r2))*log2(1-E(r2))

Inf6633=  (-0.66*log2(0.66)-0.34*log2(0.34))   +  (-0.33*log2(0.33)-0.67*log2(0.67))

Inf7525=  (-0.75*log2(0.75)-0.25*log2(0.25))   +  (-0.25*log2(0.25)-0.75*log2(0.75))

Inf9060=  (-0.90*log2(0.90)-0.10*log2(0.10))   +  (-0.60*log2(0.60)-0.40*log2(0.40))




Incertainty_Inf=[Inf6633;Inf7525;Inf9060;]
Incertainty_Inf2=[66/(33+66);75/(25+75);90/(60+90)]





