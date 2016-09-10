


function [mPdfs,l,pr,po,di,m,uniqMAP]=GirshickPredictionsWithlearntAndCardinalPrior

%strengths of the measurement density (same for likelihood), cardinal and
%learnt priors
km=[33];
kcardinal=[33];
klearnt=[33];

figure('color','w','position',[1 1 1658 457]);
colorM=[0.9 0.65 0.5;
    0.8 0 0;
    0.5 0 0];
for i=1:numel(km)
    k.m=km(i);
    k.pc=kcardinal(i);
    k.pl=klearnt(i);
    [mPdfs,l,pr,po,di,m,uniqMAP,whichplot(i)]=GirshickEs4predWithlearntAndCardinalPrior(k,colorM(i,:));
end

%annotate
legend(whichplot,'No prior','Strong cardinal prior')

