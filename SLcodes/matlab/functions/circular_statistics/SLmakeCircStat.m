%SLmakeCircStat.m
%
%       Author: Steeve Laquitaine
%         Date: 1408008
%         purpose: get descriptive stats of circular data
%
%                   - circular means, std, sem, var per condition
%                   - single conditions
%                   - number of single conditions
%         
%        usage:
%                 f0 = randsample(1:1:360,1000,'true');
%                 f1 = randsample([0.24 0.12 0.06],1000,'true');
%                 f2 = randsample([80 40 20 10],1000,'true');
%                 data = randsample(1:1:360,1000,'true');
%                 stats = SLmakeCircStat(data,f0,f1,f2)
%Description:
%    This code was previously makeSta 
%     data are in degrees.
%     f0,f1,f2 are three vectors (factors) associated with the data.

function stats = SLmakeCircStat(data,d,coh,pstd)

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
myCond = [d coh pstd];

%conditions
[conditions,subs,~] = SLuniqpair(myCond);
numConditions = size(conditions,1);
stats.conditions = conditions;
stats.conditionsSubs = subs;

%number of conditions
stats.numCondition = numConditions;

%cartesian
dataCoor = polar2cartesian(data,ones(size(data)));

%count data per condition
stats.count = nan(size(conditions,1),1);
for i = 1 : stats.numCondition 
    
    %position this condition in data set
    thiC = SLfindRow(conditions(i,:),myCond);
    
    %count
    stats.count(i) = sum(thiC);
    
    %mean, std, var, sem
    BasicStat = SLstatcircular(dataCoor(thiC==1,:));
    stats.mean(i) = BasicStat.deg.mean;    
    stats.std(i) = BasicStat.deg.std;
    stats.var(i) = BasicStat.deg.var;
    stats.sem(i) = BasicStat.deg.sem;
end




