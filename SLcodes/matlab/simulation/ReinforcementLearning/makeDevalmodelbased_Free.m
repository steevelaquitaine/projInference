
% Steeve Laquitaine
% 111028
% re-valuation by model-based/model-free controllers
% after post-training reward devaluation (e.g., pre-feeding)

numTrial=400;
% reward
r=[rand(1,numTrial/2)<0.75 rand(1,numTrial/2)<0.25]';
% action
a=ones(numTrial,1); 
% initial value
Qca=repmat(0.5,numTrial,1);
Qma=repmat(0.5,numTrial,1);
% learning rate (cached)
alpha=0.1;

for t=6:numTrial
    % values
    Qca(t)=Qca(t-1)+alpha*(r(t-1)-Qca(t-1));
    Qma(t)=r(t);
    hold all
    plot(t,Qca(t),'.-r', 'markersize',12)
    plot(t,Qma(t),'o-b')
    drawnow
end
    