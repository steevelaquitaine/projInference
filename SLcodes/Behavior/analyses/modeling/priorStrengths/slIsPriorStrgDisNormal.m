
%slIsPriorStrgDisNormal.m
%
%
% author: steeve laquitaine
%purpose: check whether prior strengths are normally distributed across 
%         subjects for the Basic Bayesian observer and the switching observer
%         p < 0.05 indicates that the null hypothesis that the data are normally 
%         distributed is rejected (a non parametric test should be used)
% 
%usage :
%
%        cd('~/Dropbox/myDropbox/Codes/projInference/data/') 
%        [~,~, data] = slcsvRead('data_PriorStrngths.csv'); 
%        [Pbo,Psw,subBOstd,subswstd,subExpstd] = slIsPriorStrgDisNormal(data);


function [Pbo,Psw,subBOstd,subswstd,subExpstd] = slIsPriorStrgDisNormal(data)

%get experimental prior strengths
subExpstd = data(1:4,:);

%get model data
subBOstd = data(5:8,:);
subswstd = data(9:12,:);

%p < 0.05 indicates that the null hypothesis that the data are normally 
%distributed is rejected (a non parametric test should be used)
for i = 1 : 4
    [~,Pbo(i)] = lillietest(subBOstd(i,:));
    [~,Psw(i)] = lillietest(subswstd(i,:));
end
  
