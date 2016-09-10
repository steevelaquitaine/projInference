

%sleyePosMeanAndCIByPrior.m
%
%  author: steeve laquitaine
%purpose : calculate mean and confidence interval of trial-averaged eye 
%          position over subjects
%
%   usage:
%
%           cd ~/.../projInference/data
%           [rowheader,colheader,data] = slcsvRead('data_meanEyePosition.csv');
%           [m,ci,condition] = sleyePosMeanAndCIByPrior(data)


function [m,ci,xy] = sleyePosMeanAndCIByPrior(data)

%sort by priors
xy80 = [data(:,1) data(:,5)];
xy40 = [data(:,2) data(:,6)];
xy20 = [data(:,3) data(:,7)];
xy10 = [data(:,4) data(:,8)];
xy = {'xy80','xy40','xy20','xy10'}';

%get mean distance and confidence interval over subjects for each prior
for i = 1 : 4    
    
    %distance of average eye position 
    %with respect to fixation point
    [~,~,~,~,~,veclen] = slgetVectorStats(eval(xy{i}));  
    
    %ci and mean
    [~,ci(i,:),m(i,:)] = slMakeCI(veclen,.95);   
    
end