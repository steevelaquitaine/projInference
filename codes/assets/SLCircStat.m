
%SLCircStat.m
%
%     Author: Steeve Laquitaine
%       Date: 131007 last modif: 140601
%
%      usage:
%        f0=randsample(1:1:360,1000,'true');
%        f1=randsample([0.24 0.12 0.06],1000,'true');
%        f2=randsample([80 40 20 10],1000,'true');
%        data=f0;
%
%        [means,stds]=SLCircStat(data,f0,f1,f2)
%
%
%        f0=randsample(5:10:355,1000,'true');
%        f1=randsample([0.24 0.12 0.06],1000,'true');
%        f2=randsample([80 40],1000,'true');
%        data=f0;
%        [means,stds,F]=SLCircStat(data,f0,f1,f2)
%
%
% Description:
%   data are in degree.
%   f0,f1,f2 are three vectors of conditions associated with the data.
%
% noteTODO:
%   add an "errorbar" option where data can be plot like with
%   errorbar in the plot for the means instead of with just dots (figure 1).

function [means,stds,F,isMeanNull] = SLCircStat(data,d,coh,pstd)


%If only one factor is input, ignore other two factors
if isempty(coh)
    coh = ones(numel(d),1);
end
if isempty(pstd)
    pstd = ones(numel(d),1);
end

%data(cartesians)
datacart = SLpolar2cartesian(data,1,'polar');

%factors 1,2,3
F.f1.i=d;
F.f1.nm='d';
F.f1.L=unique(F.f1.i);
F.f1.L=sort(F.f1.L,'ascend');
F.f1.n=numel(F.f1.L);

F.f2.i=coh;
F.f2.nm='coh';
F.f2.L=unique(F.f2.i);
F.f2.L=sort(F.f2.L,'descend');
F.f2.n=numel(F.f2.L);

F.f3.i=pstd;
F.f3.nm='pstd';
F.f3.L=unique(F.f3.i);
F.f3.L=sort(F.f3.L,'descend');
F.f3.n=numel(F.f3.L);

%positions main
for i=1:F.f1.n
    F.f1.pos(i)={find(F.f1.i==F.f1.L(i))};
end
for i=1:F.f2.n
    F.f2.pos(i)={find(F.f2.i==F.f2.L(i))};
end
for i=1:F.f3.n
    F.f3.pos(i)={find(F.f3.i==F.f3.L(i))};
end

%positions inter
for k=1:F.f1.n
    for j=1:F.f2.n
        for i=1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end

%Calculate vectorial mean and std
means=[];
stds=[];
for j=1:F.f2.n
    for k=1:F.f1.n
        for i=1:F.f3.n            
            
            stat{k,i,j} = SLcircMeanStd(datacart(F.inter.pos{k,i,j},:),'cartesian');
            means(k,i,j) = stat{k,i,j}.deg.mean;
            stds(k,i,j) = stat{k,i,j}.deg.std;
            
            
            %tag special cases when means are undefined
            %because angles cancel each other : 
            %mean coordinates are 0,0 and
            %mean angle is null
            isMeanNull(k,i,j) = stat{k,i,j}.deg.isMeanNull;
            
            %conditions
            F.f1.D3(k,i,j)=F.f1.L(k);
            F.f2.D3(k,i,j)=F.f2.L(j);
            F.f3.D3(k,i,j)=F.f3.L(i);
        end
    end
end





