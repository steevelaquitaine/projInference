
%SLbootstrapMedianStd.m
%
%author :  steeve laquitaine
%  date : 140903
%purpose: calculate bootstrapped std of data median
%
%  usage:
%
%       stdSampleMedians = SLbootstrapMedianStd(rand(10.1),1000);

function [stdSampleMedians,semSampleMedians] = SLbootstrapMedianStd(data,numboot)

%sample size
numData = numel(data);

%sample
for i = 1 : numboot
   Sample(:,i) = randsample(data,numData,'true');
end

%sample medians
SampleMedians =  nanmedian(Sample);

%sample medians std
stdSampleMedians = nanstd(SampleMedians);

%sem
SampleMedians = SLmakeColumn(SampleMedians);
semSampleMedians = sem(SampleMedians,2);