
%slMakeMeanAndStd.m
%
%
%author: steeve laquitaine
%  date: 151016
%purpose: calculate mean of mean std of means 
%         based on "the variance of the sums is the sum of the variances"
%         !!! WARNING: This is not always true. It might not be if the
%         variables are correlated !!!
%         Thus the "variance of the means is the mean of the variances"
%         converting individual means std into an std for the average of the means
%
%usage:
%
%       [meanOfmean,stdOfmeans] = slMakeMeanAndStd([1 2; 3 4],[0.1 0.2; 0.3 0.1])
%
%ref:
%http://stats.stackexchange.com/questions/21104/calculate-average-of-a-set-numbers-with-reported-standard-errors

function [meanOfmean,stdOfmeans] = slMakeMeanAndStd(means,stds)

%calculate new std based on : "The variance of the sums is the sum of the variances"
%converting individual means std into an std for the average of the means
num = size(means,1);
meanOfmean = mean(means,1);
stdOfmeans = sqrt(sum(stds.^2)/(num.^2));

