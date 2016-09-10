
      %Author: Steeve Laquitaine
        %Date: 130401
  %last modif: 130401
       %usage:
%                 f0=randsample(1:1:360,10000,'true');
%                 f1=randsample([0.24 0.12 0.06],10000,'true');
%                 f2=randsample([80 40 20 10],10000,'true');
%                 data=randsample(1:1:360,10000,'true');
%             
%                 [means,stds,F,stats] = makeStat(data',f0,f1,f2)
%Description:
    %data are in degrees.
    %f0,f1,f2 are three vectors (factors) associated with the data.

function [means,stds,F,stats] = makeStat(data,d,coh,pstd)

%check if only one factor is input, ignore other two factors
if nargin<3
    coh = ones(numel(d),1);
    pstd = ones(numel(d),1);
end

%make column vectors
data = SLmakeColumn(data);
d = SLmakeColumn(d);
coh = SLmakeColumn(coh);
pstd = SLmakeColumn(pstd);

%factors 1,2,3
F.f1.i = d;
F.f1.nm = 'd';
F.f1.L = unique(F.f1.i);
F.f1.L = sort(F.f1.L,'ascend');
F.f1.n = numel(F.f1.L);

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
for k = 1:F.f1.n
    for j = 1:F.f2.n
        for i = 1:F.f3.n
            F.inter.pos(k,i,j)=...
                {intersect( ...
                intersect(F.f1.pos{k},F.f2.pos{j}),...
                F.f3.pos{i})};
        end
    end
end
        
%mean & std
for j=1:F.f2.n
    for k=1:F.f1.n
        for i=1:F.f3.n            
            means(k,i,j) = nanmean(data(F.inter.pos{k,i,j},:));
            stds(k,i,j) = nanstd(data(F.inter.pos{k,i,j},:));
        end
    end
end



%---------------------------------
%This method is better and quicker
%---------------------------------
%I should slowly get rid of the nethod above. Too slow.
%
%descriptive stats
%conditions
[conditions,subs,~] = SLuniqpair([F.f3.i F.f2.i F.f1.i]);
numConditions = size(conditions,1);
stats.conditions = conditions;
stats.conditionsSubs = subs;
%number of conditions
stats.numCondition = numConditions;

%count data per condition
stats.count = nan(size(conditions,1),1);
for i = 1 : stats.numCondition 
    
    %position this condition in data set
    thiC = SLfindRow(conditions(i,:),[F.f3.i F.f2.i F.f1.i]);
    
    %count
    stats.count(i) = sum(thiC);
    
    %mean and std
    stats.mean(i) = nanmean(data(thiC==1));
    stats.std(i) = nanstd(data(thiC==1));
end




