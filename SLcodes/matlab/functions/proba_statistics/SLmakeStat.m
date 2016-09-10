%SLmakeStat.m
%
%       Author: Steeve Laquitaine
%         Date: 1408008
%         purpose: get descriptive stats of data:
%
%                   - means and std per condition
%                   - single conditions
%                   - number of single conditions
%         
%        usage:
%                 f0 = randsample(1:1:360,1000,'true');
%                 f1 = randsample([0.24 0.12 0.06],1000,'true');
%                 f2 = randsample([80 40 20 10],1000,'true');
%                 data = randsample(1:1:360,1000,'true');
%                 stats = SLmakeStat(data,f0,f1,f2)
%Description:
%   This code was previously makeSta 
    %data are in degrees.
    %f0,f1,f2 are three vectors (factors) associated with the data.

function stats = SLmakeStat(data,x1,x2,x3)

%check if only one factor is input, ignore other two factors

%case one factor
if nargin<3
    x2 = ones(numel(x1),1);
    x3 = ones(numel(x1),1);
end

%case two factor
if nargin<4
    x3 = ones(numel(x1),1);
end

%make column vectors
data = SLmakeColumn(data);
x1 = SLmakeColumn(x1);
x2 = SLmakeColumn(x2);
x3 = SLmakeColumn(x3);
myCond = [x1 x2 x3];

%descriptive stats
%conditions
[conditions,subs,~] = SLuniqpair(myCond);
numConditions = size(conditions,1);
stats.conditions = conditions;
stats.conditionsSubs = subs;

%number of conditions
stats.numCondition = numConditions;

%count data per condition
stats.count = nan(size(conditions,1),1);
for i = 1 : stats.numCondition 
    
    %position this condition in data set
    thiC = SLfindRow(conditions(i,:),myCond);
    
    %stat only if dataset > 1 sample
    if sum(thiC)~=1
        
        %count
        stats.count(i) = sum(thiC);
        
        %mean and std
        stats.mean(i) = nanmean(data(thiC==1));
        stats.std(i) = nanstd(data(thiC==1));
        [~,stats.confInterv{i}] = slMakeCI(data(thiC==1),0.95);
        stats.sem(i) = stats.std(i)/sqrt(stats.count(i));
    else
        %otherwise meaningless
        %count
        stats.count(i) = 1;
        stats.mean(i) = NaN;
        stats.std(i) = NaN;
        stats.confInterv{i} = [NaN NaN];
        stats.sem(i) = NaN;
    end
end





