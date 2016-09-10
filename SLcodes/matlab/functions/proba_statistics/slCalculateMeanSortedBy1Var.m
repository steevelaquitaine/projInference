%slCalculateMeanSortedBy1Var.m
%
%       Author: Steeve Laquitaine
%         Date: 1408008
%      purpose: get mean of data sorted by one variable
%
%        usage:
%                 f0 = randsample([0.24 0.12 0.06],1000,'true');
%                 data = randsample(1:1:360,1000,'true');
%                 [meanbyCond,conds] = slCalculateMeanSortedBy1Var(data,f0)
%
%        inputs :
%               data and x1 must be column vectors

function [meanbyCond,conds] = slCalculateMeanSortedBy1Var(data,x1)

%locate each condition
[conds,~,cond_location] = unique(x1);

%average data for each condition
meanbyCond = accumarray(cond_location,data,[length(conds) 1],@mean);






