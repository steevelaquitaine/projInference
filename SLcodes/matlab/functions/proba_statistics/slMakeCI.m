
%slMakeCI.m
%
%
% author: steeve laquitaine
%purpose: calculate confidence interval
%
%  usage:
%
%       [err_margin,ci,m,s] = slMakeCI(rand(40,1),.95)
%
%Description:
%
%       err_margin = z-val * s/sqrt(n)
%       ci = m +- err_margin 
%
%with s the sample std, n sample size, and the z-value (e.g., 1.96 for 95% 
%confidence level), m the sample mean
%
%
%input:
%      y  : data
%   level : 0.95 or 0.99
%
%Report statistics APA style:
%       see http://blog.apastyle.org/apastyle/2010/06/formatting-statistics-using-brackets.html

function [err_margin,ci,m,s] = slMakeCI(y,level)

y = SLmakeColumn(y);
m = nanmean(y);
s = nanstd(y);
n = length(y);

if n==1
    fprintf('%s \n','(slMakeCI) CI over one sample is meaningless. Skip this dataset.')
    err_margin = NaN
    ci = NaN
    m = NaN
    s = NaN
    return
end
    
if level==0.95
    err_margin = 1.96*s/sqrt(n);   
    ci = [m-err_margin m+err_margin];
end

if level==0.99
    err_margin = 2.58*s/sqrt(n);   
    ci = [m-err_margin m+err_margin];
end
