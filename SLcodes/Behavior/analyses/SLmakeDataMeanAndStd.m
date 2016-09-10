
% SLmakeDataMeanAndStd.m
%
%     author: steeve laquitaine
%       date: 140801
%
%    purpose: calculate data mean and std per conditions (3 conditions)
%
%      usage:
%           3 vector conditions (motion direction, coherence, std of distrib.)
%           
%           disp=repmat([1 2],3*4,1); disp=disp(:);
%           coh=repmat([1 .5],6,2); coh=coh(:);
%           pstd=repmat([10 20],3,4); pstd=pstd(:);
%           data=rand(24,1);
%           [meanData,stdData,dataCond] = SLmakeDataMeanAndStd(data,disp,...
%           coh,pstd,[],'vonMisesPrior')


%data mean and std for each experimental condition
function [meanData,stdData,dataCond] = SLmakeDataMeanAndStd(data,disp,...
    coh,pstd,priormodes,priorShape)

%(case von Mises prior)
%----------------------
%Calculate data mean and std for each experimental condition (disp,coh,pstd)
if strcmp(priorShape,'vonMisesPrior')    
    %data are sorted by experimental conditions
    [meanData,stdData,myF] = SLCircStat(data,disp,coh,pstd);
    meanData=meanData(:);
    stdData=stdData(:);
    myF1=myF.f1.D3(:);
    myF2=myF.f2.D3(:);
    myF3=myF.f3.D3(:);
    dataCond=[myF3 myF2 myF1];
end

%remove missing conditions (NaN) in the data
pos=~isnan(meanData);
meanData=meanData(pos);
stdData=stdData(pos);
dataCond=dataCond(pos,:);
