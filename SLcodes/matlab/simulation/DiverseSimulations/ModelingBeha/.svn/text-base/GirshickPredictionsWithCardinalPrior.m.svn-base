

%What we learn is that when the prior becomes too strong MAPs estimates are
%attracted such that MAPs become very similar. To be able to differentiate
%each MAPs in this situation requires a high direction resolution. Thus it 
%is slow to compute. Finding the likelihood of any given MAP is thus still
%problematic because we need to interpolate a smooth MAPs pdf and the MAP
%pdf can only be smooth if we use high resolution...

function [mPdfs,l,pr,po,di,m,uniqMAP]=GirshickPredictionsWithCardinalPrior
km=[167 16 4];
kp=[100 100 100];
figure('color','w','position',[1 1 1658 457]);
colorM=[0.9 0.65 0.5;
    0.8 0 0;
    0.5 0 0];
for i=1:numel(kp)
    k.m=km(i);
    k.p=kp(i);
    [mPdfs,l,pr,po,di,m,uniqMAP,whichplot(i)]=GirshickEs4predWithCardinalPrior(k,colorM(i,:));
end

%annotate
legend(whichplot,'No prior','Strong cardinal prior')


