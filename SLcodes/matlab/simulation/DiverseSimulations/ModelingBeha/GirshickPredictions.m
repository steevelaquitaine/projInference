

%What we learn is that when the prior becomes too strong MAPs estimates are
%attracted such that MAPs become very similar. To be able to differentiate
%each MAPs in this situation requires a high direction resolution. Thus it 
%is slow to compute. Finding the likelihood of any given MAP is thus still
%problematic because we need to interpolate a smooth MAPs pdf and the MAP
%pdf can only be smooth if we use high resolution...

function [mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs]=GirshickPredictions
km=[30 30 30];
kp=[10 30 40];
figure('color','w');
colorM=flipud(colormap('Jet'));
for i=1:numel(kp)
    k.m=km(i);
    k.p=kp(i);
    j=i*8+5;
    [mPdfs,l,pr,po,di,m,MAPfull,MAPpdfs]=GirshickEs4pred(k,colorM(j,:));
end



